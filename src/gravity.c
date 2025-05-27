#include "simple_rebound.h"

// Barnes-Hut树节点结构
typedef struct TreeNode {
    double x, y, z;           // 节点中心位置
    double size;              // 节点大小（立方体边长）
    double mass;              // 节点总质量
    double com_x, com_y, com_z; // 质心位置
    int particle_index;       // 如果是叶节点，存储粒子索引
    int is_leaf;              // 是否为叶节点
    struct TreeNode* children[8]; // 八个子节点（八叉树）
} TreeNode;


// 创建树节点
TreeNode* create_tree_node(double x, double y, double z, double size) {
    TreeNode* node = (TreeNode*)malloc(sizeof(TreeNode));
    if (!node) return NULL;
    
    node->x = x;
    node->y = y;
    node->z = z;
    node->size = size;
    node->mass = 0.0;
    node->com_x = 0.0;
    node->com_y = 0.0;
    node->com_z = 0.0;
    node->particle_index = -1;
    node->is_leaf = 1;
    
    for (int i = 0; i < 8; i++) {
        node->children[i] = NULL;
    }
    
    return node;
}

// 销毁树节点
void destroy_tree_node(TreeNode* node) {
    if (!node) return;
    
    for (int i = 0; i < 8; i++) {
        if (node->children[i]) {
            destroy_tree_node(node->children[i]);
        }
    }
    
    free(node);
}

// 确定粒子应该在哪个子节点中
int get_octant(TreeNode* node, Particle* particle) {
    int octant = 0;
    
    if (particle->x > node->x) octant |= 1;
    if (particle->y > node->y) octant |= 2;
    if (particle->z > node->z) octant |= 4;
    
    return octant;
}

// 获取子节点的中心位置
void get_child_center(TreeNode* parent, int octant, double* cx, double* cy, double* cz) {
    double half_size = parent->size * 0.25;
    
    *cx = parent->x + ((octant & 1) ? half_size : -half_size);
    *cy = parent->y + ((octant & 2) ? half_size : -half_size);
    *cz = parent->z + ((octant & 4) ? half_size : -half_size);
}

// 向树中插入粒子
void insert_particle(TreeNode* node, Simulation* sim, int particle_index) {
    Particle* particle = &sim->particles[particle_index];
    
    // 更新节点的质量和质心
    double total_mass = node->mass + particle->m;
    if (total_mass > 0) {
        node->com_x = (node->com_x * node->mass + particle->x * particle->m) / total_mass;
        node->com_y = (node->com_y * node->mass + particle->y * particle->m) / total_mass;
        node->com_z = (node->com_z * node->mass + particle->z * particle->m) / total_mass;
    }
    node->mass = total_mass;
    
    // 如果是空节点，直接插入
    if (node->is_leaf && node->particle_index == -1) {
        node->particle_index = particle_index;
        return;
    }
    
    // 如果是叶节点但已有粒子，需要分裂
    if (node->is_leaf && node->particle_index != -1) {
        int existing_particle = node->particle_index;
        node->particle_index = -1;
        node->is_leaf = 0;
        
        // 重新插入原有粒子
        int octant = get_octant(node, &sim->particles[existing_particle]);
        double cx, cy, cz;
        get_child_center(node, octant, &cx, &cy, &cz);
        
        if (!node->children[octant]) {
            node->children[octant] = create_tree_node(cx, cy, cz, node->size * 0.5);
        }
        insert_particle(node->children[octant], sim, existing_particle);
    }
    
    // 插入新粒子到适当的子节点
    int octant = get_octant(node, particle);
    double cx, cy, cz;
    get_child_center(node, octant, &cx, &cy, &cz);
    
    if (!node->children[octant]) {
        node->children[octant] = create_tree_node(cx, cy, cz, node->size * 0.5);
    }
    insert_particle(node->children[octant], sim, particle_index);
}

// 构建Barnes-Hut树
TreeNode* build_tree(Simulation* sim) {
    if (sim->N == 0) return NULL;
    
    // 找到所有粒子的边界
    double min_x = sim->particles[0].x, max_x = sim->particles[0].x;
    double min_y = sim->particles[0].y, max_y = sim->particles[0].y;
    double min_z = sim->particles[0].z, max_z = sim->particles[0].z;
    
    for (int i = 1; i < sim->N; i++) {
        if (sim->particles[i].x < min_x) min_x = sim->particles[i].x;
        if (sim->particles[i].x > max_x) max_x = sim->particles[i].x;
        if (sim->particles[i].y < min_y) min_y = sim->particles[i].y;
        if (sim->particles[i].y > max_y) max_y = sim->particles[i].y;
        if (sim->particles[i].z < min_z) min_z = sim->particles[i].z;
        if (sim->particles[i].z > max_z) max_z = sim->particles[i].z;
    }
    
    // 计算根节点的大小和中心
    double size_x = max_x - min_x;
    double size_y = max_y - min_y;
    double size_z = max_z - min_z;
    double size = fmax(fmax(size_x, size_y), size_z) * 1.1; // 稍微放大以确保包含所有粒子
    
    double center_x = (min_x + max_x) * 0.5;
    double center_y = (min_y + max_y) * 0.5;
    double center_z = (min_z + max_z) * 0.5;
    
    // 创建根节点
    TreeNode* root = create_tree_node(center_x, center_y, center_z, size);
    if (!root) return NULL;
    
    // 插入所有粒子
    for (int i = 0; i < sim->N; i++) {
        insert_particle(root, sim, i);
    }
    
    return root;
}

// 计算树节点对粒子的引力（带软化）
void calculate_force_from_node(TreeNode* node, Particle* particle, Simulation* sim, 
                              double* fx, double* fy, double* fz, double softening) {
    if (!node || node->mass == 0) return;
    
    // 计算距离
    double dx = node->com_x - particle->x;
    double dy = node->com_y - particle->y;
    double dz = node->com_z - particle->z;
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);
    
    // 避免自作用力
    if (r < 1e-10) return;
    
    // 如果是叶节点或者满足开角条件，直接计算力
    if (node->is_leaf || (node->size / r) < sim->theta) {
        // 如果是叶节点且包含当前粒子，跳过（避免自作用力）
        if (node->is_leaf && node->particle_index != -1) {
            Particle* node_particle = &sim->particles[node->particle_index];
            if (node_particle == particle) return;
        }
        
        // 应用软化
        double r2_soft = r2 + softening*softening;
        double force_magnitude = sim->G * particle->m * node->mass / r2_soft;
        double unit_x = dx / r;
        double unit_y = dy / r;
        double unit_z = dz / r;
        
        *fx += force_magnitude * unit_x;
        *fy += force_magnitude * unit_y;
        *fz += force_magnitude * unit_z;
    } else {
        // 递归计算子节点的力
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                calculate_force_from_node(node->children[i], particle, sim, fx, fy, fz, softening);
            }
        }
    }
}

// 主引力计算函数
void calculate_gravity(Simulation* sim) {
    if (!sim) return;
    
    // 根据设置的引力计算方法选择算法
    if (sim->gravity_method == GRAVITY_TREE) {
        calculate_gravity_tree(sim);
    } else {
        calculate_gravity_direct(sim);
    }
}

// 直接N²引力计算（带软化，软化系数为0时等同于不带软化）
void calculate_gravity_direct(Simulation* sim) {
    if (!sim) return;
    
    // 使用模拟中设置的软化系数
    double softening = sim->softening;
    
    // 首先清零所有粒子的加速度
    for (int i = 0; i < sim->N; i++) {
        sim->particles[i].ax = 0.0;
        sim->particles[i].ay = 0.0;
        sim->particles[i].az = 0.0;
    }
    
    // 计算每对粒子之间的引力（带软化参数）
    for (int i = 0; i < sim->N; i++) {
        for (int j = i + 1; j < sim->N; j++) {
            Particle* pi = &sim->particles[i];
            Particle* pj = &sim->particles[j];
            
            // 计算距离向量
            double dx = pj->x - pi->x;
            double dy = pj->y - pi->y;
            double dz = pj->z - pi->z;
            
            // 计算软化距离: r_soft² = r² + ε²
            double r2 = dx*dx + dy*dy + dz*dz;
            double r2_soft = r2 + softening*softening;
            double r = sqrt(r2);
            
            // 避免除零问题
            if (r < 1e-10) continue;
            
            // 计算引力大小: F = G*m1*m2/r_soft²
            double force_magnitude = sim->G * pi->m * pj->m / r2_soft;
            
            // 计算单位向量（使用原始距离）
            double unit_x = dx / r;
            double unit_y = dy / r;
            double unit_z = dz / r;
            
            // 计算加速度: a = F/m
            double ax_i = force_magnitude * unit_x / pi->m;
            double ay_i = force_magnitude * unit_y / pi->m;
            double az_i = force_magnitude * unit_z / pi->m;
            
            double ax_j = -force_magnitude * unit_x / pj->m;
            double ay_j = -force_magnitude * unit_y / pj->m;
            double az_j = -force_magnitude * unit_z / pj->m;
            
            // 累加加速度
            pi->ax += ax_i;
            pi->ay += ay_i;
            pi->az += az_i;
            
            pj->ax += ax_j;
            pj->ay += ay_j;
            pj->az += az_j;
        }
    }
}

// 树算法引力计算（带软化，软化系数为0时等同于不带软化）
void calculate_gravity_tree(Simulation* sim) {
    if (!sim || sim->N == 0) return;
    double softening = sim->softening;
    
    // 首先清零所有粒子的加速度
    for (int i = 0; i < sim->N; i++) {
        sim->particles[i].ax = 0.0;
        sim->particles[i].ay = 0.0;
        sim->particles[i].az = 0.0;
    }
    
    // 构建Barnes-Hut树
    TreeNode* root = build_tree(sim);
    if (!root) {
        printf("Gravity Method: Failed to build tree, falling back to direct method\n");
        calculate_gravity_direct(sim);
        return;
    }
    
    // 对每个粒子计算引力
    for (int i = 0; i < sim->N; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        
        // 计算来自树的总引力
        calculate_force_from_node(root, &sim->particles[i], sim, &fx, &fy, &fz, softening);
        
        // 转换为加速度
        sim->particles[i].ax = fx / sim->particles[i].m;
        sim->particles[i].ay = fy / sim->particles[i].m;
        sim->particles[i].az = fz / sim->particles[i].m;
    }
    
    // 清理树结构
    destroy_tree_node(root);
}
