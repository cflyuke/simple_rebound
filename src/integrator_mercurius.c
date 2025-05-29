#include "simple_rebound.h"
#include <math.h>

// 初始化Mercurius积分器
void mercurius_init(Simulation* sim) {
    if (!sim) return;
    
    // 设置默认参数
    sim->ri_mercurius.hill_switch_factor = 3.0;    // 希尔半径切换因子
    sim->ri_mercurius.step_switch_factor = 3.0;    // 步长相对距离切换半径
    sim->ri_mercurius.ias15_dt = 1e-12;       // IAS15时间步长（更小的步长）
    sim->ri_mercurius.current_integrator = 0;      // 默认使用WHFast
    sim->ri_mercurius.original_dt = sim->dt;       // 保存原始时间步长
    
    // 分配内存
    mercurius_alloc(sim);
    
    // 初始化子积分器
    whfast_init(sim);
    // IAS15会在第一次调用时自动初始化
}

// 重置Mercurius积分器
void mercurius_reset(Simulation* sim) {
    if (!sim) return;
    
    free(sim->ri_mercurius.distances);
    free(sim->ri_mercurius.close_encounter_flags);
    sim->ri_mercurius.distances = NULL;
    sim->ri_mercurius.close_encounter_flags = NULL;
    sim->ri_mercurius.N_allocated = 0;
    
    // 重置子积分器
    whfast_reset(sim);
    ias15_reset(sim);
}

// 分配内存
void mercurius_alloc(Simulation* sim) {
    if (!sim) return;
    
    int N = sim->N;
    int N_pairs = N * (N - 1) / 2;  // 粒子对数量
    
    if (N_pairs > sim->ri_mercurius.N_allocated) {
        sim->ri_mercurius.distances = realloc(sim->ri_mercurius.distances, sizeof(double) * N_pairs);
        sim->ri_mercurius.close_encounter_flags = realloc(sim->ri_mercurius.close_encounter_flags, sizeof(int) * N_pairs);
        sim->ri_mercurius.N_allocated = N_pairs;
        
        // 初始化标志
        for (int i = 0; i < N_pairs; i++) {
            sim->ri_mercurius.close_encounter_flags[i] = 0;
        }
    }
}

// 计算希尔半径
double mercurius_calculate_hill_radius(Particle* p1, Particle* p2, double central_mass) {
    if (!p1 || !p2) return 0.0;
    
    // 计算两个粒子到中心天体的距离
    double r1 = sqrt(p1->x * p1->x + p1->y * p1->y + p1->z * p1->z);
    double r2 = sqrt(p2->x * p2->x + p2->y * p2->y + p2->z * p2->z);
    
    // 使用较小的轨道半径
    double a = (r1 < r2) ? r1 : r2;
    
    // 计算约化质量
    double mu = (p1->m + p2->m) / central_mass;
    
    // 希尔半径公式: R_hill = a * (mu/3)^(1/3)
    return a * pow(mu / 3.0, 1.0/3.0);
}

// 检查近距离遭遇
int mercurius_check_close_encounters(Simulation* sim) {
    if (!sim || sim->N < 2) return 0;
    
    mercurius_alloc(sim);
    
    int has_close_encounter = 0;
    int pair_index = 0;
    
    // 假设第一个粒子是中心天体（质量最大）
    double central_mass = sim->particles[0].m;
    int central_index = 0;
    for (int i = 1; i < sim->N; i++) {
        if (sim->particles[i].m > central_mass) {
            central_mass = sim->particles[i].m;
            central_index = i;
        }
    }
    
    // 检查所有粒子对
    for (int i = 0; i < sim->N; i++) {
        for (int j = i + 1; j < sim->N; j++) {
            if (pair_index >= sim->ri_mercurius.N_allocated) break;
            
            Particle* p1 = &sim->particles[i];
            Particle* p2 = &sim->particles[j];
            
            // 计算粒子间距离
            double dx = p1->x - p2->x;
            double dy = p1->y - p2->y;
            double dz = p1->z - p2->z;
            double distance = sqrt(dx*dx + dy*dy + dz*dz);
            
            sim->ri_mercurius.distances[pair_index] = distance;
            
            // 如果其中一个是中心天体，跳过希尔半径检查
            if (i == central_index || j == central_index) {
                sim->ri_mercurius.close_encounter_flags[pair_index] = 0;
            } else {
                // 计算希尔半径
                double hill_radius = mercurius_calculate_hill_radius(p1, p2, central_mass);
                double switch_distance = sim->ri_mercurius.hill_switch_factor * hill_radius;

                double vx = p1->vx - p2->vx;
                double vy = p1->vy - p2->vy;
                double vz = p1->vz - p2->vz;
                double v = sqrt(vx*vx + vy*vy + vz*vz);
                double step_distance = v * sim->dt * sim->ri_mercurius.step_switch_factor;
                
                // 检查是否发生近距离遭遇
                if (distance < MIN(switch_distance, step_distance)) {
                    sim->ri_mercurius.close_encounter_flags[pair_index] = 1;
                    has_close_encounter = 1;
                } else {
                    sim->ri_mercurius.close_encounter_flags[pair_index] = 0;
                }
            }
            
            pair_index++;
        }
    }
    
    return has_close_encounter;
}

// 主Mercurius积分器函数
void integrator_mercurius(Simulation* sim) {
    if (!sim) return;
    
    // 初始化（如果需要）
    if (sim->ri_mercurius.N_allocated == 0) {
        mercurius_init(sim);
    }
    
    // 检查近距离遭遇
    int has_close_encounter = mercurius_check_close_encounters(sim);
    
    // 根据是否有近距离遭遇选择积分器
    if (has_close_encounter) {
        // 有近距离遭遇，使用IAS15积分器并减小时间步长
        if (sim->ri_mercurius.current_integrator != 1) {
            // 切换到IAS15
            sim->ri_mercurius.current_integrator = 1;
            sim->dt = sim->ri_mercurius.ias15_dt;
            
            // 确保WHFast同步（如果之前在使用）
            if (sim->ri_whfast.is_synchronized == 0) {
                whfast_synchronize(sim);
            }
        }
        
        // 使用IAS15积分器
        integrator_ias15(sim);
        
    } else {
        // 没有近距离遭遇，使用WHFast积分器
        if (sim->ri_mercurius.current_integrator != 0) {
            // 切换到WHFast
            sim->ri_mercurius.current_integrator = 0;
            sim->dt = sim->ri_mercurius.original_dt;
        }
        
        // 使用WHFast积分器
        integrator_whfast(sim);
    }
}
