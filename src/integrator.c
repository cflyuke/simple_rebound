#include "simple_rebound.h"

// 主积分器步进函数
void integrator_step(Simulation* sim) {
    if (!sim) return;
    
    switch (sim->integrator) {
        case INTEGRATOR_LEAPFROG:
            integrator_leapfrog(sim);
            apply_boundary_conditions(sim);
            sim->t += sim->dt;
            break;
        case INTEGRATOR_IAS15:
            integrator_ias15(sim);
            apply_boundary_conditions(sim);
            break;
        case INTEGRATOR_WHFAST:
            integrator_whfast(sim);
            apply_boundary_conditions(sim);
            break;
        case INTEGRATOR_MERCURIUS:
            integrator_mercurius(sim);
            apply_boundary_conditions(sim);
            break;
        default:
            integrator_ias15(sim); 
            apply_boundary_conditions(sim);
            break;
    }
}

// Leapfrog积分器（辛积分器，能量守恒较好）
void integrator_leapfrog(Simulation* sim) {
    if (!sim) return;
    
    // 第一步：更新速度半步
    calculate_gravity(sim);
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        p->vx += p->ax * sim->dt * 0.5;
        p->vy += p->ay * sim->dt * 0.5;
        p->vz += p->az * sim->dt * 0.5;
    }
    
    // 第二步：更新位置全步
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        p->x += p->vx * sim->dt;
        p->y += p->vy * sim->dt;
        p->z += p->vz * sim->dt;
    }
    
    // 第三步：重新计算加速度并更新速度另外半步
    calculate_gravity(sim);
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        p->vx += p->ax * sim->dt * 0.5;
        p->vy += p->ay * sim->dt * 0.5;
        p->vz += p->az * sim->dt * 0.5;
    }
}
