#include "simple_rebound.h"

// 创建新的模拟
Simulation* sim_create(void) {
    Simulation* sim = (Simulation*)malloc(sizeof(Simulation));
    if (!sim) {
        fprintf(stderr, "Error: Failed to allocate memory for simulation\n");
        return NULL;
    }
    
    sim->N = 0;
    sim->N_allocated = 100;  // 初始分配100个粒子的空间
    sim->particles = (Particle*)malloc(sim->N_allocated * sizeof(Particle));
    if (!sim->particles) {
        fprintf(stderr, "Error: Failed to allocate memory for particles\n");
        free(sim);
        return NULL;
    }
    
    sim->t = 0.0;
    sim->dt = 0.001;
    sim->exact_finish_time = 0;  // 默认不使用精确时间控制
    sim->G = G_DEFAULT;
    sim->softening = 0.0;           // 默认软化系数为0
    sim->theta = 0.5;               // 默认Barnes-Hut开角参数为0.5
    sim->integrator = INTEGRATOR_IAS15;
    sim->collision = COLLISION_MERGE;
    sim->collision_detection = COLLISION_DETECTION_DIRECT;  // 默认使用直接检测方法
    sim->gravity_method = GRAVITY_DIRECT;  // 默认使用直接方法
    sim->boundary = BOUNDARY_OPEN;  // 默认开放边界
    sim->boundary_size = 100.0;     // 默认边界大小
    sim->energy = 0.0;
    sim->momentum[0] = sim->momentum[1] = sim->momentum[2] = 0.0;
    
    // 初始化IAS15积分器
    memset(&sim->ri_ias15, 0, sizeof(IAS15_integrator));
    sim->ri_ias15.epsilon = 1e-9;
    sim->ri_ias15.min_dt = 0;
    sim->ri_ias15.adaptive_mode = 2;
    sim->ri_ias15.iterations_max_exceeded = 0;
    
    // 初始化WHFast积分器
    memset(&sim->ri_whfast, 0, sizeof(WHFast_integrator));
    sim->ri_whfast.corrector = 0;  // 默认无校正器
    sim->ri_whfast.coordinates = WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;  // 使用democratic heliocentric坐标系
    sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
    sim->ri_whfast.safe_mode = 0;  // 默认DKD方案
    sim->ri_whfast.keep_unsynchronized = 0;
    sim->ri_whfast.is_synchronized = 1;
    sim->ri_whfast.N_allocated = 0;
    sim->ri_whfast.timestep_warning = 0;
    sim->ri_whfast.recalculate_coordinates_but_not_synchronized_warning = 0;
    sim->ri_whfast.p_jh = NULL;
    
    // 初始化Mercurius积分器
    memset(&sim->ri_mercurius, 0, sizeof(Mercurius_integrator));
    sim->ri_mercurius.hill_switch_factor = 3.0;
    sim->ri_mercurius.step_switch_factor = 3.0;
    sim->ri_mercurius.ias15_dt = 1e-12;
    sim->ri_mercurius.current_integrator = 0;
    sim->ri_mercurius.original_dt = sim->dt;
    sim->ri_mercurius.distances = NULL;
    sim->ri_mercurius.close_encounter_flags = NULL;
    sim->ri_mercurius.N_allocated = 0;
    
    return sim;
}

// 销毁模拟
void sim_destroy(Simulation* sim) {
    if (sim) {
        if (sim->particles) {
            free(sim->particles);
        }
        // 清理IAS15积分器
        ias15_reset(sim);
        
        // 清理WHFast积分器
        if (sim->ri_whfast.p_jh) {
            free(sim->ri_whfast.p_jh);
            sim->ri_whfast.p_jh = NULL;
        }
        sim->ri_whfast.N_allocated = 0;
        
        // 清理Mercurius积分器
        mercurius_reset(sim);
        
        free(sim);
    }
}

// 添加粒子
void sim_add_particle(Simulation* sim, double m, double x, double y, double z, 
                     double vx, double vy, double vz, double r) {
    if (!sim) return;
    
    // 检查是否需要扩展数组
    if (sim->N >= sim->N_allocated) {
        sim->N_allocated *= 2;
        sim->particles = (Particle*)realloc(sim->particles, 
                                          sim->N_allocated * sizeof(Particle));
        if (!sim->particles) {
            fprintf(stderr, "Error: Failed to reallocate memory for particles\n");
            return;
        }
    }
    
    Particle* p = &sim->particles[sim->N];
    p->x = x; p->y = y; p->z = z;
    p->vx = vx; p->vy = vy; p->vz = vz;
    p->ax = 0.0; p->ay = 0.0; p->az = 0.0;
    p->m = m;
    p->r = r;
    p->id = sim->N;
    
    sim->N++;
    
    // 标记需要重新计算坐标
    sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
}

// 移除粒子
void sim_remove_particle(Simulation* sim, int index) {
    if (!sim || index < 0 || index >= sim->N) return;
    
    // 将最后一个粒子移动到被删除粒子的位置
    if (index < sim->N - 1) {
        sim->particles[index] = sim->particles[sim->N - 1];
        sim->particles[index].id = index;  // 更新ID
    }
    sim->N--;
    
    // 标记需要重新计算坐标
    sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
}

// 设置软化系数
void sim_set_softening(Simulation* sim, double softening) {
    if (!sim) return;
    sim->softening = softening;
}

// 设置Barnes-Hut开角参数
void sim_set_theta(Simulation* sim, double theta) {
    if (!sim) return;
    sim->theta = theta;
}

// 设置引力计算方法
void sim_set_gravity_method(Simulation* sim, GravityType method) {
    if (!sim) return;
    sim->gravity_method = method;
}

// 设置碰撞检测方法
void sim_set_collision_detection(Simulation* sim, CollisionDetectionType method) {
    if (!sim) return;
    sim->collision_detection = method;
}

// 设置边界条件
void sim_set_boundary(Simulation* sim, BoundaryType boundary, double size) {
    if (!sim) return;
    sim->boundary = boundary;
    sim->boundary_size = size;
}

// 设置精确时间控制
void sim_set_exact_finish_time(Simulation* sim, int exact) {
    if (!sim) return;
    sim->exact_finish_time = exact;
}

// 应用边界条件
void apply_boundary_conditions(Simulation* sim) {
    if (!sim) return;
    
    if (sim->boundary == BOUNDARY_OPEN) {
        // 开放边界：删除离开边界的粒子
        double half_size = sim->boundary_size * 0.5;
        
        for (int i = sim->N - 1; i >= 0; i--) {  // 从后往前遍历，避免删除时索引问题
            Particle* p = &sim->particles[i];
            
            // 检查粒子是否超出边界
            if (fabs(p->x) > half_size || fabs(p->y) > half_size || fabs(p->z) > half_size) {
                sim_remove_particle(sim, i);
            }
        }
    } else {
        // 其他边界条件
        for (int i = 0; i < sim->N; i++) {
            Particle* p = &sim->particles[i];
            
            switch (sim->boundary) {
                case BOUNDARY_PERIODIC:
                    apply_periodic_boundary(sim, p);
                    break;
                case BOUNDARY_REFLECTIVE:
                    apply_reflective_boundary(sim, p);
                    break;
                default:
                    break;
            }
        }
    }
}

// 应用周期性边界条件
void apply_periodic_boundary(Simulation* sim, Particle* p) {
    if (!sim || !p) return;
    
    double half_size = sim->boundary_size * 0.5;
    
    // X方向
    if (p->x > half_size) {
        p->x -= sim->boundary_size;
    } else if (p->x < -half_size) {
        p->x += sim->boundary_size;
    }
    
    // Y方向
    if (p->y > half_size) {
        p->y -= sim->boundary_size;
    } else if (p->y < -half_size) {
        p->y += sim->boundary_size;
    }
    
    // Z方向
    if (p->z > half_size) {
        p->z -= sim->boundary_size;
    } else if (p->z < -half_size) {
        p->z += sim->boundary_size;
    }
}

// 应用反射边界条件
void apply_reflective_boundary(Simulation* sim, Particle* p) {
    if (!sim || !p) return;
    
    double half_size = sim->boundary_size * 0.5;
    
    // X方向
    if (p->x > half_size) {
        p->x = 2.0 * half_size - p->x;
        p->vx = -p->vx;
    } else if (p->x < -half_size) {
        p->x = -2.0 * half_size - p->x;
        p->vx = -p->vx;
    }
    
    // Y方向
    if (p->y > half_size) {
        p->y = 2.0 * half_size - p->y;
        p->vy = -p->vy;
    } else if (p->y < -half_size) {
        p->y = -2.0 * half_size - p->y;
        p->vy = -p->vy;
    }
    
    // Z方向
    if (p->z > half_size) {
        p->z = 2.0 * half_size - p->z;
        p->vz = -p->vz;
    } else if (p->z < -half_size) {
        p->z = -2.0 * half_size - p->z;
        p->vz = -p->vz;
    }
}

// 计算距离
double distance(Particle* p1, Particle* p2) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// 计算系统总能量
double calculate_energy(Simulation* sim) {
    if (!sim) return 0.0;
    
    double kinetic = 0.0;
    double potential = 0.0;
    
    // 计算动能
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        kinetic += 0.5 * p->m * (p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
    }
    
    // 计算势能
    for (int i = 0; i < sim->N; i++) {
        for (int j = i + 1; j < sim->N; j++) {
            Particle* pi = &sim->particles[i];
            Particle* pj = &sim->particles[j];
            double r = distance(pi, pj);
            if (r > 0) {
                // 考虑软化系数
                if (sim->softening > 0) {
                    double r_soft = sqrt(r*r + sim->softening*sim->softening);
                    potential -= sim->G * pi->m * pj->m / r_soft;
                } else {
                    potential -= sim->G * pi->m * pj->m / r;
                }
            }
        }
    }
    
    sim->energy = kinetic + potential;
    return sim->energy;
}

// 计算系统总动量
void calculate_momentum(Simulation* sim) {
    if (!sim) return;
    
    sim->momentum[0] = sim->momentum[1] = sim->momentum[2] = 0.0;
    
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        sim->momentum[0] += p->m * p->vx;
        sim->momentum[1] += p->m * p->vy;
        sim->momentum[2] += p->m * p->vz;
    }
}

// 打印模拟信息
void print_simulation_info(Simulation* sim) {
    if (!sim) return;
    
    printf("Simulation Info:\n");
    printf("  Time: %.6e\n", sim->t);
    printf("  Particles: %d\n", sim->N);
    printf("  Time step: %.6e\n", sim->dt);
    printf("  G: %.6e\n", sim->G);
    printf("  Softening: %.6e\n", sim->softening);
    printf("  Theta: %.6e\n", sim->theta);
    
    const char* integrator_str;
    switch (sim->integrator) {
        case INTEGRATOR_LEAPFROG:
            integrator_str = "Leapfrog";
            break;
        case INTEGRATOR_IAS15:
            integrator_str = "IAS15";
            break;
        case INTEGRATOR_WHFAST:
            integrator_str = "WHFast";
            break;
        case INTEGRATOR_MERCURIUS:
            integrator_str = "Mercurius";
            break;
        default:
            integrator_str = "Unknown";
            break;
    }
    printf("  Integrator: %s\n", integrator_str);
    
    printf("  Gravity method: %s\n", 
           sim->gravity_method == GRAVITY_TREE ? "Tree (Barnes-Hut)" : "Direct");
    printf("  Collision detection: %s\n", 
           sim->collision_detection == COLLISION_DETECTION_SPATIAL ? "Spatial" : "Direct");
    
    const char* boundary_str;
    switch (sim->boundary) {
        case BOUNDARY_OPEN:
            boundary_str = "Open";
            break;
        case BOUNDARY_PERIODIC:
            boundary_str = "Periodic";
            break;
        case BOUNDARY_REFLECTIVE:
            boundary_str = "Reflective";
            break;
        default:
            boundary_str = "Unknown";
            break;
    }
    printf("  Boundary: %s", boundary_str);
    if (sim->boundary != BOUNDARY_OPEN) {
        printf(" (size: %.6e)", sim->boundary_size);
    }
    printf("\n");
    
    printf("  Energy: %.6e\n", calculate_energy(sim));
    calculate_momentum(sim);
    printf("  Momentum: (%.6e, %.6e, %.6e)\n", 
           sim->momentum[0], sim->momentum[1], sim->momentum[2]);
}
