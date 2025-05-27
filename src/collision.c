#include "simple_rebound.h"

// 主碰撞检测函数
void check_collisions(Simulation* sim) {
    if (!sim || sim->collision == COLLISION_NONE) return;
    
    // 根据设置的碰撞检测方法选择算法
    if (sim->collision_detection == COLLISION_DETECTION_SPATIAL) {
        check_collisions_spatial(sim);
    } else {
        check_collisions_direct(sim);
    }
}

// 直接碰撞检测（O(N²)方法）
void check_collisions_direct(Simulation* sim) {
    if (!sim || sim->collision == COLLISION_NONE) return;
    
    // 检查所有粒子对是否发生碰撞
    for (int i = 0; i < sim->N; i++) {
        for (int j = i + 1; j < sim->N; j++) {
            Particle* pi = &sim->particles[i];
            Particle* pj = &sim->particles[j];
            
            // 计算粒子间距离
            double dist = distance(pi, pj);
            
            // 检查是否发生碰撞（距离小于半径之和）
            if (dist < (pi->r + pj->r)) {
                resolve_collision(sim, i, j);
                // 碰撞后粒子数量可能改变，需要调整循环
                if (sim->collision == COLLISION_MERGE && j < sim->N) {
                    j--; // 重新检查当前位置的粒子
                }
            }
        }
    }
}

// 碰撞处理函数
void resolve_collision(Simulation* sim, int i, int j) {
    if (!sim || i < 0 || j < 0 || i >= sim->N || j >= sim->N || i == j) return;
    sim->ri_whfast.recalculate_coordinates_this_timestep = 1; 
    switch (sim->collision) {
        case COLLISION_MERGE:
            merge_particles(sim, i, j);
            break;
        case COLLISION_BOUNCE:
            bounce_particles(sim, i, j);
            break;
        default:
            break;
    }
}

// 粒子合并（完全非弹性碰撞）
void merge_particles(Simulation* sim, int i, int j) {
    if (!sim || i < 0 || j < 0 || i >= sim->N || j >= sim->N || i == j) return;
    
    Particle* pi = &sim->particles[i];
    Particle* pj = &sim->particles[j];
    
    // 保存总质量和总动量
    double total_mass = pi->m + pj->m;
    double total_px = pi->m * pi->vx + pj->m * pj->vx;
    double total_py = pi->m * pi->vy + pj->m * pj->vy;
    double total_pz = pi->m * pi->vz + pj->m * pj->vz;
    
    // 计算质心位置（按质量加权）
    double cm_x = (pi->m * pi->x + pj->m * pj->x) / total_mass;
    double cm_y = (pi->m * pi->y + pj->m * pj->y) / total_mass;
    double cm_z = (pi->m * pi->z + pj->m * pj->z) / total_mass;
    
    // 计算新的速度（动量守恒）
    double new_vx = total_px / total_mass;
    double new_vy = total_py / total_mass;
    double new_vz = total_pz / total_mass;
    
    // 计算新的半径（假设密度相同）
    double volume_i = (4.0/3.0) * M_PI * pi->r * pi->r * pi->r;
    double volume_j = (4.0/3.0) * M_PI * pj->r * pj->r * pj->r;
    double new_volume = volume_i + volume_j;
    double new_radius = pow(new_volume * 3.0 / (4.0 * M_PI), 1.0/3.0);
    
    // 更新第一个粒子的属性
    pi->x = cm_x;
    pi->y = cm_y;
    pi->z = cm_z;
    pi->vx = new_vx;
    pi->vy = new_vy;
    pi->vz = new_vz;
    pi->m = total_mass;
    pi->r = new_radius;
    
    // 移除第二个粒子
    sim_remove_particle(sim, j);
}

// 粒子弹性碰撞
void bounce_particles(Simulation* sim, int i, int j) {
    if (!sim || i < 0 || j < 0 || i >= sim->N || j >= sim->N || i == j) return;
    
    Particle* pi = &sim->particles[i];
    Particle* pj = &sim->particles[j];
    
    // 计算相对位置和速度
    double dx = pj->x - pi->x;
    double dy = pj->y - pi->y;
    double dz = pj->z - pi->z;
    
    double dvx = pj->vx - pi->vx;
    double dvy = pj->vy - pi->vy;
    double dvz = pj->vz - pi->vz;
    
    // 计算距离
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    if (dist < 1e-10) return; // 避免除零
    
    // 单位法向量
    double nx = dx / dist;
    double ny = dy / dist;
    double nz = dz / dist;
    
    // 相对速度在法向上的分量
    double dvn = dvx*nx + dvy*ny + dvz*nz;
    
    // 如果粒子正在分离，不需要处理碰撞
    if (dvn > 0) return;
    
    // 计算冲量（弹性碰撞）
    double impulse = 2.0 * dvn / (pi->m + pj->m);
    
    // 更新速度
    pi->vx += impulse * pj->m * nx;
    pi->vy += impulse * pj->m * ny;
    pi->vz += impulse * pj->m * nz;
    
    pj->vx -= impulse * pi->m * nx;
    pj->vy -= impulse * pi->m * ny;
    pj->vz -= impulse * pi->m * nz;
    
    // 分离粒子以避免重叠
    double overlap = (pi->r + pj->r) - dist;
    if (overlap > 0) {
        double separation = overlap * 0.5 + 1e-10; // 稍微多分离一点
        
        pi->x -= separation * nx;
        pi->y -= separation * ny;
        pi->z -= separation * nz;
        
        pj->x += separation * nx;
        pj->y += separation * ny;
        pj->z += separation * nz;
    }
}

// 网格单元结构
typedef struct GridCell {
    int* particle_indices;
    int count;
    int capacity;
} GridCell;

// 空间网格结构
typedef struct SpatialGrid {
    GridCell* cells;
    int grid_size_x, grid_size_y, grid_size_z;
    double cell_size;
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
} SpatialGrid;

// 创建空间网格
static SpatialGrid* create_spatial_grid(Simulation* sim) {
    if (!sim || sim->N == 0) return NULL;
    
    SpatialGrid* grid = (SpatialGrid*)malloc(sizeof(SpatialGrid));
    if (!grid) return NULL;
    
    // 找到粒子的边界框
    grid->min_x = grid->max_x = sim->particles[0].x;
    grid->min_y = grid->max_y = sim->particles[0].y;
    grid->min_z = grid->max_z = sim->particles[0].z;
    
    double max_radius = sim->particles[0].r;
    
    for (int i = 1; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        if (p->x < grid->min_x) grid->min_x = p->x;
        if (p->x > grid->max_x) grid->max_x = p->x;
        if (p->y < grid->min_y) grid->min_y = p->y;
        if (p->y > grid->max_y) grid->max_y = p->y;
        if (p->z < grid->min_z) grid->min_z = p->z;
        if (p->z > grid->max_z) grid->max_z = p->z;
        if (p->r > max_radius) max_radius = p->r;
    }
    
    // 扩展边界以包含粒子半径
    grid->min_x -= max_radius;
    grid->min_y -= max_radius;
    grid->min_z -= max_radius;
    grid->max_x += max_radius;
    grid->max_y += max_radius;
    grid->max_z += max_radius;
    
    // 设置网格单元大小（至少是最大粒子直径的2倍）
    grid->cell_size = fmax(2.0 * max_radius, 
                          fmax((grid->max_x - grid->min_x) / 50.0,
                               fmax((grid->max_y - grid->min_y) / 50.0,
                                    (grid->max_z - grid->min_z) / 50.0)));
    
    // 计算网格维度
    grid->grid_size_x = (int)ceil((grid->max_x - grid->min_x) / grid->cell_size) + 1;
    grid->grid_size_y = (int)ceil((grid->max_y - grid->min_y) / grid->cell_size) + 1;
    grid->grid_size_z = (int)ceil((grid->max_z - grid->min_z) / grid->cell_size) + 1;
    
    // 限制网格大小以避免内存过度使用
    if (grid->grid_size_x > 100) grid->grid_size_x = 100;
    if (grid->grid_size_y > 100) grid->grid_size_y = 100;
    if (grid->grid_size_z > 100) grid->grid_size_z = 100;
    
    // 分配网格单元
    int total_cells = grid->grid_size_x * grid->grid_size_y * grid->grid_size_z;
    grid->cells = (GridCell*)calloc(total_cells, sizeof(GridCell));
    if (!grid->cells) {
        free(grid);
        return NULL;
    }
    
    // 初始化每个单元
    for (int i = 0; i < total_cells; i++) {
        grid->cells[i].capacity = 10; // 初始容量
        grid->cells[i].particle_indices = (int*)malloc(grid->cells[i].capacity * sizeof(int));
        grid->cells[i].count = 0;
        if (!grid->cells[i].particle_indices) {
            // 清理已分配的内存
            for (int j = 0; j < i; j++) {
                free(grid->cells[j].particle_indices);
            }
            free(grid->cells);
            free(grid);
            return NULL;
        }
    }
    
    return grid;
}

// 销毁空间网格
static void destroy_spatial_grid(SpatialGrid* grid) {
    if (!grid) return;
    
    int total_cells = grid->grid_size_x * grid->grid_size_y * grid->grid_size_z;
    for (int i = 0; i < total_cells; i++) {
        free(grid->cells[i].particle_indices);
    }
    free(grid->cells);
    free(grid);
}

// 获取粒子所在的网格单元索引
static int get_grid_index(SpatialGrid* grid, double x, double y, double z) {
    int gx = (int)((x - grid->min_x) / grid->cell_size);
    int gy = (int)((y - grid->min_y) / grid->cell_size);
    int gz = (int)((z - grid->min_z) / grid->cell_size);
    
    // 边界检查
    if (gx < 0) gx = 0;
    if (gy < 0) gy = 0;
    if (gz < 0) gz = 0;
    if (gx >= grid->grid_size_x) gx = grid->grid_size_x - 1;
    if (gy >= grid->grid_size_y) gy = grid->grid_size_y - 1;
    if (gz >= grid->grid_size_z) gz = grid->grid_size_z - 1;
    
    return gz * grid->grid_size_x * grid->grid_size_y + gy * grid->grid_size_x + gx;
}

// 向网格单元添加粒子
static void add_particle_to_cell(GridCell* cell, int particle_index) {
    if (cell->count >= cell->capacity) {
        cell->capacity *= 2;
        cell->particle_indices = (int*)realloc(cell->particle_indices, 
                                              cell->capacity * sizeof(int));
        if (!cell->particle_indices) return; // 内存分配失败
    }
    
    cell->particle_indices[cell->count] = particle_index;
    cell->count++;
}

// 获取相邻网格单元的索引
static void get_neighbor_cells(SpatialGrid* grid, int gx, int gy, int gz, 
                              int* neighbor_indices, int* neighbor_count) {
    *neighbor_count = 0;
    
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                int nx = gx + dx;
                int ny = gy + dy;
                int nz = gz + dz;
                
                if (nx >= 0 && nx < grid->grid_size_x &&
                    ny >= 0 && ny < grid->grid_size_y &&
                    nz >= 0 && nz < grid->grid_size_z) {
                    
                    int index = nz * grid->grid_size_x * grid->grid_size_y + 
                               ny * grid->grid_size_x + nx;
                    neighbor_indices[(*neighbor_count)++] = index;
                }
            }
        }
    }
}

// 高效碰撞检测（空间分割方法）
void check_collisions_spatial(Simulation* sim) {
    if (!sim || sim->collision == COLLISION_NONE || sim->N < 2) return;
    
    // 对于少量粒子，直接使用O(N²)方法更高效
    if (sim->N < 50) {
        check_collisions(sim);
        return;
    }
    
    // 创建空间网格
    SpatialGrid* grid = create_spatial_grid(sim);
    if (!grid) {
        // 如果网格创建失败，回退到简单方法
        printf("Collision Detection: failed to create grid, back to the direct method.");
        check_collisions(sim);
        return;
    }
    
    // 将粒子分配到网格单元
    for (int i = 0; i < sim->N; i++) {
        Particle* p = &sim->particles[i];
        int cell_index = get_grid_index(grid, p->x, p->y, p->z);
        add_particle_to_cell(&grid->cells[cell_index], i);
    }
    
    // 检查每个粒子与其相邻单元中的粒子的碰撞
    for (int i = 0; i < sim->N; i++) {
        Particle* pi = &sim->particles[i];
        
        // 获取粒子所在的网格坐标
        int gx = (int)((pi->x - grid->min_x) / grid->cell_size);
        int gy = (int)((pi->y - grid->min_y) / grid->cell_size);
        int gz = (int)((pi->z - grid->min_z) / grid->cell_size);
        
        // 边界检查
        if (gx < 0) gx = 0;
        if (gy < 0) gy = 0;
        if (gz < 0) gz = 0;
        if (gx >= grid->grid_size_x) gx = grid->grid_size_x - 1;
        if (gy >= grid->grid_size_y) gy = grid->grid_size_y - 1;
        if (gz >= grid->grid_size_z) gz = grid->grid_size_z - 1;
        
        // 获取相邻单元
        int neighbor_indices[27]; // 最多27个相邻单元（3x3x3）
        int neighbor_count;
        get_neighbor_cells(grid, gx, gy, gz, neighbor_indices, &neighbor_count);
        
        // 检查与相邻单元中粒子的碰撞
        for (int n = 0; n < neighbor_count; n++) {
            GridCell* cell = &grid->cells[neighbor_indices[n]];
            
            for (int c = 0; c < cell->count; c++) {
                int j = cell->particle_indices[c];
                
                // 避免重复检查和自碰撞
                if (j <= i) continue;
                
                // 检查粒子是否仍然存在（可能在之前的碰撞中被移除）
                if (i >= sim->N || j >= sim->N) continue;
                
                Particle* pj = &sim->particles[j];
                
                // 计算粒子间距离
                double dist = distance(pi, pj);
                
                // 检查是否发生碰撞
                if (dist < (pi->r + pj->r)) {
                    resolve_collision(sim, i, j);
                    
                    // 碰撞后粒子数量可能改变，需要重新构建网格
                    if (sim->collision == COLLISION_MERGE) {
                        destroy_spatial_grid(grid);
                        // 递归调用以重新构建网格
                        check_collisions_spatial(sim);
                        return;
                    }
                }
            }
        }
    }
    
    // 清理网格
    destroy_spatial_grid(grid);
}

// 碰撞预测（连续碰撞检测）
double predict_collision_time(Particle* p1, Particle* p2) {
    // 计算两个粒子的相对位置和速度
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    double dz = p2->z - p1->z;
    
    double dvx = p2->vx - p1->vx;
    double dvy = p2->vy - p1->vy;
    double dvz = p2->vz - p1->vz;
    
    // 求解二次方程: |r(t)|² = (r_sum)²
    // 其中 r(t) = r0 + v*t, r_sum = r1 + r2
    double a = dvx*dvx + dvy*dvy + dvz*dvz;
    double b = 2.0 * (dx*dvx + dy*dvy + dz*dvz);
    double c = dx*dx + dy*dy + dz*dz - (p1->r + p2->r) * (p1->r + p2->r);
    
    // 判别式
    double discriminant = b*b - 4*a*c;
    
    if (discriminant < 0 || a == 0) {
        return -1.0; // 无碰撞
    }
    
    // 求解时间（取较小的正值）
    double t1 = (-b - sqrt(discriminant)) / (2*a);
    double t2 = (-b + sqrt(discriminant)) / (2*a);
    
    if (t1 > 0) return t1;
    if (t2 > 0) return t2;
    
    return -1.0; // 无未来碰撞
}
