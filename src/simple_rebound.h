#ifndef SIMPLE_REBOUND_H
#define SIMPLE_REBOUND_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 常量定义
#define G_DEFAULT 6.67430e-11
#define MAX_PARTICLES 10000

#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define MIN(a, b) ((a) > (b) ? (b) : (a))

// WHFast相关常量
#define WHFAST_NMAX_QUART 64    // 四次求解器最大迭代次数
#define WHFAST_NMAX_NEWT  32    // Newton方法最大迭代次数


// 粒子结构
typedef struct {
    double x, y, z;     // 位置
    double vx, vy, vz;  // 速度
    double ax, ay, az;  // 加速度
    double m;           // 质量
    double r;           // 半径
    int id;             // 粒子ID
}Particle;

// 积分器类型
typedef enum {
    INTEGRATOR_LEAPFROG,
    INTEGRATOR_IAS15,
    INTEGRATOR_WHFAST,
    INTEGRATOR_MERCURIUS
} IntegratorType;

// 碰撞处理类型
typedef enum {
    COLLISION_NONE,
    COLLISION_MERGE,
    COLLISION_BOUNCE
} CollisionType;

// 引力计算方法类型
typedef enum {
    GRAVITY_DIRECT,
    GRAVITY_TREE
} GravityType;

// 碰撞检测方法类型
typedef enum {
    COLLISION_DETECTION_DIRECT,
    COLLISION_DETECTION_SPATIAL
} CollisionDetectionType;

// 边界类型
typedef enum {
    BOUNDARY_OPEN,      // 开放边界（无边界）
    BOUNDARY_PERIODIC,  // 周期性边界
    BOUNDARY_REFLECTIVE // 反射边界
} BoundaryType;

// IAS15积分器相关结构
typedef struct {
    double* p0;
    double* p1;
    double* p2;
    double* p3;
    double* p4;
    double* p5;
    double* p6;
} IAS15_dp7;

typedef struct {
    IAS15_dp7 g;        // g coefficients
    IAS15_dp7 b;        // b coefficients  
    IAS15_dp7 csb;      // compensated summation for b
    IAS15_dp7 e;        // e coefficients
    IAS15_dp7 br;       // backup b coefficients
    IAS15_dp7 er;       // backup e coefficients
    double* at;         // temporary acceleration storage
    double* x0;         // initial positions
    double* v0;         // initial velocities
    double* a0;         // initial accelerations
    double* csx;        // compensated summation for x
    double* csv;        // compensated summation for v
    double* csa0;       // compensated summation for a0
    int* map;           // particle mapping
    int N_allocated;    // allocated size
    int N_allocated_map; // allocated map size
    double epsilon;     // tolerance for adaptive timestep
    double min_dt;      // minimum timestep
    double dt_last_done; // last completed timestep
    int adaptive_mode;  // adaptive timestep mode
    int iterations_max_exceeded; // counter for max iterations exceeded
} IAS15_integrator;

// WHFast坐标系类型
typedef enum {
    WHFAST_COORDINATES_JACOBI = 0,                      // Jacobi坐标系（默认）
    WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC = 1,      // Democratic Heliocentric坐标系
    WHFAST_COORDINATES_WHDS = 2,                        // WHDS坐标系
    WHFAST_COORDINATES_BARYCENTRIC = 3,                 // 重心坐标系
} WHFastCoordinates;

// WHFast积分器结构
typedef struct {
    unsigned int corrector;                             // 辛校正器阶数：0（默认-无校正器），3，5，7，11，17
    WHFastCoordinates coordinates;                      // 坐标系类型
    unsigned int recalculate_coordinates_this_timestep; // 1：从惯性坐标重新计算坐标
    unsigned int safe_mode;                             // 0：不启用安全模式（默认）；1：启用安全模式
    unsigned int keep_unsynchronized;                   // 1：同步后继续非同步状态
    
    // 内部使用
    Particle* p_jh;                                     // Jacobi/heliocentric/WHDS坐标
    unsigned int is_synchronized;                       // 同步状态标志
    unsigned int N_allocated;                           // 已分配大小
    unsigned int timestep_warning;                      // 时间步长警告计数器
    unsigned int recalculate_coordinates_but_not_synchronized_warning; // 坐标重计算但未同步警告
} WHFast_integrator;

// Mercurius积分器结构
typedef struct {
    double hill_switch_factor;      // 希尔半径切换因子，默认为3.0
    double whfast_dt_factor;        // WHFast时间步长因子，默认为1.0
    double ias15_dt_factor;         // IAS15时间步长因子，默认为0.1
    double* distances;              // 粒子间距离缓存
    int* close_encounter_flags;     // 近距离遭遇标志
    int N_allocated;                // 已分配大小
    int current_integrator;         // 当前使用的积分器 (0=WHFast, 1=IAS15)
    double original_dt;             // 原始时间步长
} Mercurius_integrator;

// 模拟结构
typedef struct {
    Particle* particles;
    int N;                      // 粒子数量
    int N_allocated;            // 分配的粒子数量
    double t;                   // 当前时间
    double dt;                  // 时间步长
    int exact_finish_time;      // 是否需要精确到达结束时间
    double G;                   // 引力常数
    double softening;           // 软化系数，默认为0
    double theta;               // Barnes-Hut开角参数，默认为0.5
    IntegratorType integrator;  // 积分器类型
    CollisionType collision;    // 碰撞处理类型
    CollisionDetectionType collision_detection; // 碰撞检测方法
    GravityType gravity_method; // 引力计算方法
    BoundaryType boundary;      // 边界类型
    double boundary_size;       // 边界大小（立方体边长）
    double energy;              // 系统能量
    double momentum[3];         // 系统动量
    IAS15_integrator ri_ias15;  // IAS15积分器数据
    WHFast_integrator ri_whfast; // WHFast积分器数据
    Mercurius_integrator ri_mercurius; // Mercurius积分器数据
} Simulation;

// 函数声明

// 模拟管理
Simulation* sim_create(void);
void sim_destroy(Simulation* sim);
void sim_add_particle(Simulation* sim, double m, double x, double y, double z, 
                     double vx, double vy, double vz, double r);
void sim_remove_particle(Simulation* sim, int index);
void sim_set_softening(Simulation* sim, double softening);
void sim_set_theta(Simulation* sim, double theta);
void sim_set_gravity_method(Simulation* sim, GravityType method);
void sim_set_collision_detection(Simulation* sim, CollisionDetectionType method);
void sim_set_boundary(Simulation* sim, BoundaryType boundary, double size);
void sim_set_exact_finish_time(Simulation* sim, int exact);

// 积分器
void integrator_step(Simulation* sim);
void integrator_leapfrog(Simulation* sim);
void integrator_ias15(Simulation* sim);
void integrator_whfast(Simulation* sim);
void integrator_mercurius(Simulation* sim);

// IAS15积分器相关函数
void ias15_alloc(Simulation* sim);
void ias15_reset(Simulation* sim);
int ias15_step(Simulation* sim);

// WHFast积分器相关函数
void whfast_init(Simulation* sim);
void whfast_reset(Simulation* sim);
void whfast_kepler_solver(Simulation* sim, Particle* p_j, double M, int i, double _dt);
void whfast_from_inertial(Simulation* sim);
void whfast_to_inertial(Simulation* sim);
void whfast_interaction_step(Simulation* sim, double _dt);
void whfast_jump_step(Simulation* sim, double _dt);
void whfast_kepler_step(Simulation* sim, double _dt);
void whfast_com_step(Simulation* sim, double _dt);
void whfast_synchronize(Simulation* sim);

// Mercurius积分器相关函数
void mercurius_init(Simulation* sim);
void mercurius_reset(Simulation* sim);
void mercurius_alloc(Simulation* sim);
int mercurius_check_close_encounters(Simulation* sim);
double mercurius_calculate_hill_radius(Particle* p1, Particle* p2, double central_mass);

// 坐标变换函数
void particles_transform_inertial_to_democraticheliocentric_posvel(const Particle* const particles, Particle* const p_h, const unsigned int N, const unsigned int N_active);
void particles_transform_democraticheliocentric_to_inertial_posvel(Particle* const particles, const Particle* const p_h, const unsigned int N, const unsigned int N_active);
void particles_transform_democraticheliocentric_to_inertial_pos(Particle* const particles, const Particle* const p_h, const unsigned int N, const unsigned int N_active);

// 引力计算
void calculate_gravity(Simulation* sim);
void calculate_gravity_direct(Simulation* sim);
void calculate_gravity_tree(Simulation* sim);

// 碰撞检测
void check_collisions(Simulation* sim);
void check_collisions_direct(Simulation* sim);
void check_collisions_spatial(Simulation* sim);
void resolve_collision(Simulation* sim, int i, int j);
void merge_particles(Simulation* sim, int i, int j);
void bounce_particles(Simulation* sim, int i, int j);
double predict_collision_time(Particle* p1, Particle* p2);

// 边界处理
void apply_boundary_conditions(Simulation* sim);
void apply_periodic_boundary(Simulation* sim, Particle* p);
void apply_reflective_boundary(Simulation* sim, Particle* p);

// 粒子工具函数
Particle create_particle(double m, double x, double y, double z, 
                        double vx, double vy, double vz, double r);
Particle copy_particle(const Particle* src);
double particle_kinetic_energy(const Particle* p);
void particle_momentum(const Particle* p, double* px, double* py, double* pz);
void particle_angular_momentum(const Particle* p, double* lx, double* ly, double* lz);
double particle_speed(const Particle* p);
double particle_distance_from_origin(const Particle* p);
void set_particle_orbit(Particle* p, double central_mass, double semi_major_axis, 
                       double eccentricity, double inclination, double longitude_of_node,
                       double argument_of_periapsis, double mean_anomaly, double G);
double relative_velocity(const Particle* p1, const Particle* p2);
double escape_velocity(const Particle* p1, const Particle* p2, double G);
void print_particle_info(const Particle* p);
int is_particle_valid(const Particle* p);

// 工具函数
double calculate_energy(Simulation* sim);
void calculate_momentum(Simulation* sim);
double distance(Particle* p1, Particle* p2);
void print_simulation_info(Simulation* sim);

#endif // SIMPLE_REBOUND_H
