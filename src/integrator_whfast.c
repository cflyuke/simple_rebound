/**
 * @file    integrator_whfast.c
 * @brief   WHFast Kepler solver implementation for simple_rebound
 * @details This file implements the core Kepler solver from WHFast integrator
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simple_rebound.h"

// 快速绝对值函数
static inline double fastabs(double x) {
    return (x > 0.) ? x : -x;
}

// 快速阶乘倒数查找表
static const double invfactorial[35] = {
    1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 
    1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 
    1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 
    1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 
    1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 
    1./15511210043330985984000000., 1./403291461126605635584000000., 
    1./10888869450418352160768000000., 1./304888344611713860501504000000., 
    1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 
    1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 
    1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.
};

/**
 * 计算Stumpff函数 c0, c1, c2, c3, c4, c5
 * 用于处理不同类型的轨道（椭圆、抛物线、双曲线）
 */
static void stumpff_cs(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fastabs(z) > 0.1) {
        z = z / 4.;
        n++;
    }
    const int nmax = 15;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np = nmax-2; np >= 5; np -= 2) {
        c_odd  = invfactorial[np]    - z * c_odd;
        c_even = invfactorial[np-1]  - z * c_even;
    }
    cs[5] = c_odd;
    cs[4] = c_even;
    cs[3] = invfactorial[3]  - z * cs[5];
    cs[2] = invfactorial[2]  - z * cs[4];
    cs[1] = invfactorial[1]  - z * cs[3];
    for (; n > 0; n--) { 
        z = z * 4.;
        cs[5] = (cs[5] + cs[4] + cs[3] * cs[2]) * 0.0625;
        cs[4] = (1. + cs[1]) * cs[3] * 0.125;
        cs[3] = 1./6. - z * cs[5];
        cs[2] = 0.5 - z * cs[4];
        cs[1] = 1. - z * cs[3];
    }
    cs[0] = invfactorial[0]  - z * cs[2];
}

/**
 * 计算简化的Stumpff函数 c0, c1, c2, c3
 */
static void stumpff_cs3(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fabs(z) > 0.1) {
        z = z / 4.;
        n++;
    }
    const int nmax = 13;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np = nmax-2; np >= 3; np -= 2) {
        c_odd  = invfactorial[np]    - z * c_odd;
        c_even = invfactorial[np-1]  - z * c_even;
    }
    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1]  - z * c_odd;
    cs[0] = invfactorial[0]  - z * c_even;
    for (; n > 0; n--) { 
        cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
        cs[2] = cs[1] * cs[1] * 0.5;
        cs[1] = cs[0] * cs[1];
        cs[0] = 2. * cs[0] * cs[0] - 1.;
    }
}

/**
 * 计算Stiefel G函数（完整版本）
 */
static void stiefel_Gs(double *restrict Gs, double beta, double X) {
    double X2 = X * X;
    stumpff_cs(Gs, beta * X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    double _pow = X2 * X;
    Gs[3] *= _pow; 
    _pow *= X;
    Gs[4] *= _pow; 
    _pow *= X;
    Gs[5] *= _pow; 
    return;
}

/**
 * 计算Stiefel G函数（简化版本，只计算G0-G3）
 */
static void stiefel_Gs3(double *restrict Gs, double beta, double X) {
    double X2 = X * X;
    stumpff_cs3(Gs, beta * X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    Gs[3] *= X2 * X;
    return;
}

/**
 * WHFast Kepler求解器 - 核心函数
 * 求解单个粒子在中心引力场中的Keplerian运动
 * 
 * @param sim 模拟结构指针
 * @param p_j 粒子数组（Jacobi坐标系）
 * @param M 中心天体的引力参数 (G*m_central)
 * @param i 粒子索引
 * @param _dt 时间步长
 */
void whfast_kepler_solver(Simulation* sim, Particle* p_j, double M, int i, double _dt) {
    const Particle p1 = p_j[i];

    // 计算初始轨道参数
    const double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
    const double r0i = 1./r0;
    const double v2 = p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
    const double beta = 2.*M*r0i - v2;  // 轨道能量参数
    const double eta0 = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;  // 径向速度
    const double zeta0 = M - beta*r0;
    
    double X;
    double Gs[6]; 
    double invperiod = 0;  // 只用于椭圆轨道
    double X_per_period = NAN; // 只用于椭圆轨道，NAN触发Newton方法用于双曲线轨道

    if (beta > 0.) {
        // 椭圆轨道
        const double sqrt_beta = sqrt(beta);
        invperiod = sqrt_beta*beta/(2.*M_PI*M);
        X_per_period = 2.*M_PI/sqrt_beta;
        
        // 检查时间步长警告
        if (fabs(_dt)*invperiod > 1. && sim->ri_whfast.timestep_warning == 0) {
            sim->ri_whfast.timestep_warning++;
            printf("Warning: WHFast convergence issue. Timestep is larger than at least one orbital period.\n");
        }
        
        // 初始猜测（二阶）
        const double dtr0i = _dt*r0i;
        X = dtr0i * (1. - dtr0i*eta0*0.5*r0i);
    } else {
        // 双曲线轨道
        X = 0.; // 初始猜测
    }

    unsigned int converged = 0;
    double oldX = X; 

    // 执行一步Newton迭代
    stiefel_Gs3(Gs, beta, X);
    const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
    double ri = 1./(r0 + eta0Gs1zeta0Gs2);
    X = ri*(X*eta0Gs1zeta0Gs2 - eta0*Gs[2] - zeta0*Gs[3] + _dt);

    // 根据估计的步长选择求解器
    // 注意：对于双曲线轨道，这使用Newton方法
    if(fastabs(X-oldX) > 0.01*X_per_period) {
        // 四次求解器
        // 线性初始猜测
        X = beta*_dt/M;
        double prevX[WHFAST_NMAX_QUART+1];
        for(int n_lag = 1; n_lag < WHFAST_NMAX_QUART; n_lag++) {
            stiefel_Gs3(Gs, beta, X);
            const double f = r0*X + eta0*Gs[2] + zeta0*Gs[3] - _dt;
            const double fp = r0 + eta0*Gs[1] + zeta0*Gs[2];
            const double fpp = eta0*Gs[0] + zeta0*Gs[1];
            const double denom = fp + sqrt(fabs(16.*fp*fp - 20.*f*fpp));
            X = (X*denom - 5.*f)/denom;
            for(int j = 1; j < n_lag; j++) {
                if(X == prevX[j]) {
                    // 收敛，退出
                    n_lag = WHFAST_NMAX_QUART;
                    converged = 1;
                    break;
                }
            }
            prevX[n_lag] = X;
        }
        const double eta0Gs1zeta0Gs2_final = eta0*Gs[1] + zeta0*Gs[2];
        ri = 1./(r0 + eta0Gs1zeta0Gs2_final);
    } else {
        // Newton方法
        double oldX2 = NAN;             
        for (int n_hg = 1; n_hg < WHFAST_NMAX_NEWT; n_hg++) {
            oldX2 = oldX;
            oldX = X;
            stiefel_Gs3(Gs, beta, X);
            const double eta0Gs1zeta0Gs2_iter = eta0*Gs[1] + zeta0*Gs[2];
            ri = 1./(r0 + eta0Gs1zeta0Gs2_iter);
            X = ri*(X*eta0Gs1zeta0Gs2_iter - eta0*Gs[2] - zeta0*Gs[3] + _dt);

            if (X == oldX || X == oldX2) {
                // 收敛，退出
                converged = 1;
                break; 
            }
        }
    }

    // 如果求解器不工作，回退到二分法
    if (converged == 0) { 
        double X_min, X_max;
        if (beta > 0.) {
            // 椭圆轨道
            X_min = X_per_period * floor(_dt*invperiod);
            X_max = X_min + X_per_period;
        } else {
            // 双曲线轨道
            double h2 = r0*r0*v2 - eta0*eta0;
            double q = h2/M/(1. + sqrt(1. - h2*beta/(M*M)));
            double vq = copysign(sqrt(h2)/q, _dt);
            // X_max和X_min对应dt/r_min和dt/r_max
            // 这些是在此时间步长内可达到的
            // r_max = vq*_dt+r0
            // r_min = 近心点
            X_min = _dt/(fastabs(vq*_dt) + r0); 
            X_max = _dt/q;
            if (_dt < 0.) {
                double temp = X_min;
                X_min = X_max;
                X_max = temp;
            }
        }
        X = (X_max + X_min)/2.;
        do {
            stiefel_Gs3(Gs, beta, X);
            double s = r0*X + eta0*Gs[2] + zeta0*Gs[3] - _dt;
            if (s >= 0.) {
                X_max = X;
            } else {
                X_min = X;
            }
            X = (X_max + X_min)/2.;
        } while (fastabs((X_max-X_min)) > fastabs((X_max+X_min)*1e-15));
        const double eta0Gs1zeta0Gs2_bisect = eta0*Gs[1] + zeta0*Gs[2];
        ri = 1./(r0 + eta0Gs1zeta0Gs2_bisect);
    }
    
    if (isnan(ri)) {
        // 双曲线情况下（几乎）直线运动的异常处理
        ri = 0.;
        Gs[1] = 0.;
        Gs[2] = 0.;
        Gs[3] = 0.;
    }

    // 注意：这些不是传统的f和g函数
    double f = -M*Gs[2]*r0i;
    double g = _dt - M*Gs[3];
    double fd = -M*Gs[1]*r0i*ri; 
    double gd = -M*Gs[2]*ri; 

    // 更新位置和速度
    p_j[i].x += f*p1.x + g*p1.vx;
    p_j[i].y += f*p1.y + g*p1.vy;
    p_j[i].z += f*p1.z + g*p1.vz;

    p_j[i].vx += fd*p1.x + gd*p1.vx;
    p_j[i].vy += fd*p1.y + gd*p1.vy;
    p_j[i].vz += fd*p1.z + gd*p1.vz;
}

/**
 * 初始化WHFast积分器
 */
void whfast_init(Simulation* sim) {
    if (!sim) return;
    
    // 初始化WHFast积分器参数
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
    
    // 分配内存
    if (sim->N > (int)sim->ri_whfast.N_allocated) {
        sim->ri_whfast.p_jh = realloc(sim->ri_whfast.p_jh, sim->N * sizeof(Particle));
        sim->ri_whfast.N_allocated = sim->N;
    }
}

/**
 * 重置WHFast积分器
 */
void whfast_reset(Simulation* sim) {
    if (!sim) return;
    
    sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
    sim->ri_whfast.is_synchronized = 1;
    sim->ri_whfast.timestep_warning = 0;
    sim->ri_whfast.recalculate_coordinates_but_not_synchronized_warning = 0;
}

/**
 * 从惯性坐标系转换到WHFast坐标系
 */
void whfast_from_inertial(Simulation* sim) {
    if (!sim) return;
    
    // 确保内存已分配
    if (sim->N > (int)sim->ri_whfast.N_allocated) {
        sim->ri_whfast.p_jh = realloc(sim->ri_whfast.p_jh, sim->N * sizeof(Particle));
        sim->ri_whfast.N_allocated = sim->N;
    }
    
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            particles_transform_inertial_to_democraticheliocentric_posvel(
                sim->particles, sim->ri_whfast.p_jh, sim->N, sim->N);
            break;
        default:
            // 默认复制惯性坐标
            memcpy(sim->ri_whfast.p_jh, sim->particles, sim->N * sizeof(Particle));
            break;
    }
    
    sim->ri_whfast.recalculate_coordinates_this_timestep = 0;
}

/**
 * 从WHFast坐标系转换回惯性坐标系
 */
void whfast_to_inertial(Simulation* sim) {
    if (!sim) return;
    
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            particles_transform_democraticheliocentric_to_inertial_posvel(
                sim->particles, sim->ri_whfast.p_jh, sim->N, sim->N);
            break;
        default:
            // 默认复制坐标
            memcpy(sim->particles, sim->ri_whfast.p_jh, sim->N * sizeof(Particle));
            break;
    }
}

/**
 * WHFast相互作用步骤（Kick步骤）
 */
void whfast_interaction_step(Simulation* sim, double _dt) {
    if (!sim) return;
    
    // 计算加速度
    calculate_gravity(sim);
    
    // 更新速度
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            // 对于democratic heliocentric坐标系，只更新非中心粒子的速度
            for (int i = 1; i < sim->N; i++) {
                sim->ri_whfast.p_jh[i].vx += _dt * sim->particles[i].ax;
                sim->ri_whfast.p_jh[i].vy += _dt * sim->particles[i].ay;
                sim->ri_whfast.p_jh[i].vz += _dt * sim->particles[i].az;
            }
            break;
        default:
            for (int i = 0; i < sim->N; i++) {
                sim->ri_whfast.p_jh[i].vx += _dt * sim->particles[i].ax;
                sim->ri_whfast.p_jh[i].vy += _dt * sim->particles[i].ay;
                sim->ri_whfast.p_jh[i].vz += _dt * sim->particles[i].az;
            }
            break;
    }
}

/**
 * WHFast Jump步骤（处理质心运动）
 */
void whfast_jump_step(Simulation* sim, double _dt) {
    if (!sim) return;
    
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            {
                if (sim->N < 2) return;
                
                // 计算总动量
                double px = 0, py = 0, pz = 0;
                double m0 = sim->particles[0].m;  // 中心天体质量
                
                for (int i = 1; i < sim->N; i++) {
                    double m = sim->particles[i].m;
                    px += m * sim->ri_whfast.p_jh[i].vx;
                    py += m * sim->ri_whfast.p_jh[i].vy;
                    pz += m * sim->ri_whfast.p_jh[i].vz;
                }
                
                // 更新所有粒子位置
                for (int i = 1; i < sim->N; i++) {
                    sim->ri_whfast.p_jh[i].x += _dt * (px / m0);
                    sim->ri_whfast.p_jh[i].y += _dt * (py / m0);
                    sim->ri_whfast.p_jh[i].z += _dt * (pz / m0);
                }
            }
            break;
        default:
            // 其他坐标系暂不实现
            break;
    }
}

/**
 * WHFast Kepler步骤（Drift步骤）
 */
void whfast_kepler_step(Simulation* sim, double _dt) {
    if (!sim) return;
    
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            {
                if (sim->N < 2) return;
                
                double M = sim->G * sim->particles[0].m;  // 中心天体的引力参数
                
                // 对每个非中心粒子求解Kepler问题
                for (int i = 1; i < sim->N; i++) {
                    whfast_kepler_solver(sim, sim->ri_whfast.p_jh, M, i, _dt);
                }
            }
            break;
        default:
            // 其他坐标系暂不实现
            break;
    }
}

/**
 * WHFast 质心步骤
 */
void whfast_com_step(Simulation* sim, double _dt) {
    if (!sim) return;
    
    switch (sim->ri_whfast.coordinates) {
        case WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            // 质心粒子的位置更新
            if (sim->N >= 1) {
                sim->ri_whfast.p_jh[0].x += sim->ri_whfast.p_jh[0].vx * _dt;
                sim->ri_whfast.p_jh[0].y += sim->ri_whfast.p_jh[0].vy * _dt;
                sim->ri_whfast.p_jh[0].z += sim->ri_whfast.p_jh[0].vz * _dt;
            }
            break;
        default:
            // 其他坐标系暂不实现
            break;
    }
}

/**
 * WHFast同步函数
 */
void whfast_synchronize(Simulation* sim) {
    if (!sim) return;
    
    if (sim->ri_whfast.is_synchronized == 0) {
        // 如果处于非同步状态，需要完成剩余的半步
        whfast_from_inertial(sim);
        whfast_kepler_step(sim, sim->dt / 2.0);
        whfast_com_step(sim, sim->dt / 2.0);
        whfast_to_inertial(sim);
        sim->ri_whfast.is_synchronized = 1;
    }
}

/**
 * 主WHFast积分器函数
 */
void integrator_whfast(Simulation* sim) {
    if (!sim) return;
    
    // 初始化（如果需要）
    if (sim->ri_whfast.N_allocated == 0) {
        whfast_init(sim);
    }
    
    // 如果需要重新计算坐标或处于同步状态
    if (sim->ri_whfast.recalculate_coordinates_this_timestep || sim->ri_whfast.safe_mode) {
        whfast_from_inertial(sim);
        //当前的实现不需要改变synchronized状态，因为没有step中途访问数据的需求。
    }
    
    // DKD积分方案
    // 第一个半步Drift
    whfast_kepler_step(sim, sim->dt / 2.0);
    whfast_com_step(sim, sim->dt / 2.0);
    
    // Jump步骤
    whfast_jump_step(sim, sim->dt / 2.0);
    
    // 转换回惯性坐标系进行相互作用计算
    whfast_to_inertial(sim);
    
    // Kick步骤（相互作用）
    whfast_interaction_step(sim, sim->dt);
    
    // 转换回WHFast坐标系
    whfast_from_inertial(sim);
    
    // 第二个Jump步骤
    whfast_jump_step(sim, sim->dt / 2.0);
    
    // 第二个半步Drift
    whfast_kepler_step(sim, sim->dt / 2.0);
    whfast_com_step(sim, sim->dt / 2.0);
    
    // 转换回惯性坐标系
    whfast_to_inertial(sim);
    
    // 更新时间
    sim->t += sim->dt;
}
