#include "simple_rebound.h"

// 粒子相关的工具函数

// 创建新粒子
Particle create_particle(double m, double x, double y, double z, 
                        double vx, double vy, double vz, double r) {
    Particle p;
    p.x = x; p.y = y; p.z = z;
    p.vx = vx; p.vy = vy; p.vz = vz;
    p.ax = 0.0; p.ay = 0.0; p.az = 0.0;
    p.m = m;
    p.r = r;
    p.id = -1;  // 将在添加到模拟时设置
    return p;
}

// 复制粒子
Particle copy_particle(const Particle* src) {
    Particle p;
    if (src) {
        p = *src;
    } else {
        p = create_particle(0, 0, 0, 0, 0, 0, 0, 0);
    }
    return p;
}

// 计算粒子的动能
double particle_kinetic_energy(const Particle* p) {
    if (!p) return 0.0;
    return 0.5 * p->m * (p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
}

// 计算粒子的动量
void particle_momentum(const Particle* p, double* px, double* py, double* pz) {
    if (!p) {
        *px = *py = *pz = 0.0;
        return;
    }
    *px = p->m * p->vx;
    *py = p->m * p->vy;
    *pz = p->m * p->vz;
}

// 计算粒子的角动量（相对于原点）
void particle_angular_momentum(const Particle* p, double* lx, double* ly, double* lz) {
    if (!p) {
        *lx = *ly = *lz = 0.0;
        return;
    }
    
    double px = p->m * p->vx;
    double py = p->m * p->vy;
    double pz = p->m * p->vz;
    
    *lx = p->y * pz - p->z * py;
    *ly = p->z * px - p->x * pz;
    *lz = p->x * py - p->y * px;
}

// 计算粒子速度的大小
double particle_speed(const Particle* p) {
    if (!p) return 0.0;
    return sqrt(p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
}

// 计算粒子位置的大小（距离原点）
double particle_distance_from_origin(const Particle* p) {
    if (!p) return 0.0;
    return sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
}

// 设置粒子的轨道参数（简化的开普勒轨道）
void set_particle_orbit(Particle* p, double central_mass, double semi_major_axis, 
                       double eccentricity, double inclination, double longitude_of_node,
                       double argument_of_periapsis, double mean_anomaly, double G) {
    if (!p) return;
    
    // 简化实现：假设圆轨道（eccentricity = 0）
    double r = semi_major_axis;
    double v = sqrt(G * central_mass / r);
    
    // 设置位置（在xy平面上）
    p->x = r * cos(mean_anomaly);
    p->y = r * sin(mean_anomaly);
    p->z = 0.0;
    
    // 设置速度（垂直于位置向量）
    p->vx = -v * sin(mean_anomaly);
    p->vy = v * cos(mean_anomaly);
    p->vz = 0.0;
    
    // 应用倾角旋转（简化）
    if (inclination != 0.0) {
        double cos_i = cos(inclination);
        double sin_i = sin(inclination);
        
        double new_y = p->y * cos_i - p->z * sin_i;
        double new_z = p->y * sin_i + p->z * cos_i;
        p->y = new_y;
        p->z = new_z;
        
        double new_vy = p->vy * cos_i - p->vz * sin_i;
        double new_vz = p->vy * sin_i + p->vz * cos_i;
        p->vy = new_vy;
        p->vz = new_vz;
    }
}

// 计算两个粒子的相对速度
double relative_velocity(const Particle* p1, const Particle* p2) {
    if (!p1 || !p2) return 0.0;
    
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;
    double dvz = p2->vz - p1->vz;
    
    return sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
}

// 计算粒子的逃逸速度（相对于另一个粒子）
double escape_velocity(const Particle* p1, const Particle* p2, double G) {
    if (!p1 || !p2) return 0.0;
    
    double r = distance((Particle*)p1, (Particle*)p2);
    if (r <= 0) return 0.0;
    
    return sqrt(2.0 * G * (p1->m + p2->m) / r);
}

// 打印粒子信息
void print_particle_info(const Particle* p) {
    if (!p) {
        printf("Particle: NULL\n");
        return;
    }
    
    printf("Particle ID: %d\n", p->id);
    printf("  Position: (%.6e, %.6e, %.6e)\n", p->x, p->y, p->z);
    printf("  Velocity: (%.6e, %.6e, %.6e)\n", p->vx, p->vy, p->vz);
    printf("  Acceleration: (%.6e, %.6e, %.6e)\n", p->ax, p->ay, p->az);
    printf("  Mass: %.6e\n", p->m);
    printf("  Radius: %.6e\n", p->r);
    printf("  Speed: %.6e\n", particle_speed(p));
    printf("  Kinetic Energy: %.6e\n", particle_kinetic_energy(p));
}

// 检查粒子是否有效（非NaN值）
int is_particle_valid(const Particle* p) {
    if (!p) return 0;
    
    return !isnan(p->x) && !isnan(p->y) && !isnan(p->z) &&
           !isnan(p->vx) && !isnan(p->vy) && !isnan(p->vz) &&
           !isnan(p->ax) && !isnan(p->ay) && !isnan(p->az) &&
           !isnan(p->m) && !isnan(p->r) &&
           !isinf(p->x) && !isinf(p->y) && !isinf(p->z) &&
           !isinf(p->vx) && !isinf(p->vy) && !isinf(p->vz) &&
           !isinf(p->ax) && !isinf(p->ay) && !isinf(p->az) &&
           !isinf(p->m) && !isinf(p->r);
}
