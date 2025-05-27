/**
 * @file    transformations.c
 * @brief   Coordinate transformation functions for WHFast integrator
 * @details This file implements coordinate transformations between inertial and
 *          democratic heliocentric coordinate systems used by WHFast
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simple_rebound.h"

/**
 * 从惯性坐标系转换到Democratic Heliocentric坐标系（位置和速度）
 * 
 * Democratic Heliocentric坐标系的定义：
 * - 第0个粒子（中心天体）：相对于系统质心的位置和速度
 * - 第i个粒子（i>0）：相对于中心天体的位置和速度
 * 
 * 这种坐标系的优点是：
 * 1. 保持了系统的总动量守恒
 * 2. 将多体问题分解为多个二体问题
 * 3. 适合用辛积分器求解
 * 
 * @param particles 输入：惯性坐标系中的粒子数组
 * @param p_h 输出：Democratic Heliocentric坐标系中的粒子数组
 * @param N 粒子总数
 * @param N_active 活跃粒子数（通常等于N）
 */
void particles_transform_inertial_to_democraticheliocentric_posvel(
    const Particle* const particles, 
    Particle* const p_h, 
    const unsigned int N, 
    const unsigned int N_active) {
    
    if (!particles || !p_h || N == 0) return;
    
    // 计算系统质心位置和速度
    double com_x = 0.0, com_y = 0.0, com_z = 0.0;
    double com_vx = 0.0, com_vy = 0.0, com_vz = 0.0;
    double total_mass = 0.0;
    
    for (unsigned int i = 0; i < N_active; i++) {
        double m = particles[i].m;
        total_mass += m;
        com_x += m * particles[i].x;
        com_y += m * particles[i].y;
        com_z += m * particles[i].z;
        com_vx += m * particles[i].vx;
        com_vy += m * particles[i].vy;
        com_vz += m * particles[i].vz;
    }
    
    if (total_mass > 0.0) {
        com_x /= total_mass;
        com_y /= total_mass;
        com_z /= total_mass;
        com_vx /= total_mass;
        com_vy /= total_mass;
        com_vz /= total_mass;
    }
    
    // 转换第0个粒子（中心天体）：相对于质心
    p_h[0].x = particles[0].x - com_x;
    p_h[0].y = particles[0].y - com_y;
    p_h[0].z = particles[0].z - com_z;
    p_h[0].vx = particles[0].vx - com_vx;
    p_h[0].vy = particles[0].vy - com_vy;
    p_h[0].vz = particles[0].vz - com_vz;
    p_h[0].m = particles[0].m;
    p_h[0].r = particles[0].r;
    p_h[0].id = particles[0].id;
    p_h[0].ax = particles[0].ax;
    p_h[0].ay = particles[0].ay;
    p_h[0].az = particles[0].az;
    
    // 转换其他粒子：相对于中心天体
    for (unsigned int i = 1; i < N; i++) {
        p_h[i].x = particles[i].x - particles[0].x;
        p_h[i].y = particles[i].y - particles[0].y;
        p_h[i].z = particles[i].z - particles[0].z;
        p_h[i].vx = particles[i].vx - particles[0].vx;
        p_h[i].vy = particles[i].vy - particles[0].vy;
        p_h[i].vz = particles[i].vz - particles[0].vz;
        p_h[i].m = particles[i].m;
        p_h[i].r = particles[i].r;
        p_h[i].id = particles[i].id;
        p_h[i].ax = particles[i].ax;
        p_h[i].ay = particles[i].ay;
        p_h[i].az = particles[i].az;
    }
}

/**
 * 从Democratic Heliocentric坐标系转换到惯性坐标系（位置和速度）
 * 
 * 这是上述变换的逆变换
 * 
 * @param particles 输出：惯性坐标系中的粒子数组
 * @param p_h 输入：Democratic Heliocentric坐标系中的粒子数组
 * @param N 粒子总数
 * @param N_active 活跃粒子数（通常等于N）
 */
void particles_transform_democraticheliocentric_to_inertial_posvel(
    Particle* const particles, 
    const Particle* const p_h, 
    const unsigned int N, 
    const unsigned int N_active) {
    
    if (!particles || !p_h || N == 0) return;
    
    // 计算Democratic Heliocentric坐标系中的质心
    double com_x = 0.0, com_y = 0.0, com_z = 0.0;
    double com_vx = 0.0, com_vy = 0.0, com_vz = 0.0;
    double total_mass = 0.0;
    
    // 第0个粒子的贡献
    double m0 = p_h[0].m;
    total_mass += m0;
    com_x += m0 * p_h[0].x;
    com_y += m0 * p_h[0].y;
    com_z += m0 * p_h[0].z;
    com_vx += m0 * p_h[0].vx;
    com_vy += m0 * p_h[0].vy;
    com_vz += m0 * p_h[0].vz;
    
    // 其他粒子的贡献（注意：它们的坐标是相对于第0个粒子的）
    for (unsigned int i = 1; i < N_active; i++) {
        double m = p_h[i].m;
        total_mass += m;
        com_x += m * (p_h[i].x + p_h[0].x);
        com_y += m * (p_h[i].y + p_h[0].y);
        com_z += m * (p_h[i].z + p_h[0].z);
        com_vx += m * (p_h[i].vx + p_h[0].vx);
        com_vy += m * (p_h[i].vy + p_h[0].vy);
        com_vz += m * (p_h[i].vz + p_h[0].vz);
    }
    
    if (total_mass > 0.0) {
        com_x /= total_mass;
        com_y /= total_mass;
        com_z /= total_mass;
        com_vx /= total_mass;
        com_vy /= total_mass;
        com_vz /= total_mass;
    }
    
    // 转换第0个粒子：从相对质心坐标到惯性坐标
    particles[0].x = p_h[0].x + com_x;
    particles[0].y = p_h[0].y + com_y;
    particles[0].z = p_h[0].z + com_z;
    particles[0].vx = p_h[0].vx + com_vx;
    particles[0].vy = p_h[0].vy + com_vy;
    particles[0].vz = p_h[0].vz + com_vz;
    particles[0].m = p_h[0].m;
    particles[0].r = p_h[0].r;
    particles[0].id = p_h[0].id;
    particles[0].ax = p_h[0].ax;
    particles[0].ay = p_h[0].ay;
    particles[0].az = p_h[0].az;
    
    // 转换其他粒子：从相对中心天体坐标到惯性坐标
    for (unsigned int i = 1; i < N; i++) {
        particles[i].x = p_h[i].x + particles[0].x;
        particles[i].y = p_h[i].y + particles[0].y;
        particles[i].z = p_h[i].z + particles[0].z;
        particles[i].vx = p_h[i].vx + particles[0].vx;
        particles[i].vy = p_h[i].vy + particles[0].vy;
        particles[i].vz = p_h[i].vz + particles[0].vz;
        particles[i].m = p_h[i].m;
        particles[i].r = p_h[i].r;
        particles[i].id = p_h[i].id;
        particles[i].ax = p_h[i].ax;
        particles[i].ay = p_h[i].ay;
        particles[i].az = p_h[i].az;
    }
}

/**
 * 从Democratic Heliocentric坐标系转换到惯性坐标系（仅位置）
 * 
 * 这个函数只转换位置，不转换速度，用于某些特殊情况
 * 
 * @param particles 输出：惯性坐标系中的粒子数组
 * @param p_h 输入：Democratic Heliocentric坐标系中的粒子数组
 * @param N 粒子总数
 * @param N_active 活跃粒子数（通常等于N）
 */
void particles_transform_democraticheliocentric_to_inertial_pos(
    Particle* const particles, 
    const Particle* const p_h, 
    const unsigned int N, 
    const unsigned int N_active) {
    
    if (!particles || !p_h || N == 0) return;
    
    // 计算Democratic Heliocentric坐标系中的质心位置
    double com_x = 0.0, com_y = 0.0, com_z = 0.0;
    double total_mass = 0.0;
    
    // 第0个粒子的贡献
    double m0 = p_h[0].m;
    total_mass += m0;
    com_x += m0 * p_h[0].x;
    com_y += m0 * p_h[0].y;
    com_z += m0 * p_h[0].z;
    
    // 其他粒子的贡献
    for (unsigned int i = 1; i < N_active; i++) {
        double m = p_h[i].m;
        total_mass += m;
        com_x += m * (p_h[i].x + p_h[0].x);
        com_y += m * (p_h[i].y + p_h[0].y);
        com_z += m * (p_h[i].z + p_h[0].z);
    }
    
    if (total_mass > 0.0) {
        com_x /= total_mass;
        com_y /= total_mass;
        com_z /= total_mass;
    }
    
    // 转换第0个粒子位置
    particles[0].x = p_h[0].x + com_x;
    particles[0].y = p_h[0].y + com_y;
    particles[0].z = p_h[0].z + com_z;
    
    // 转换其他粒子位置
    for (unsigned int i = 1; i < N; i++) {
        particles[i].x = p_h[i].x + particles[0].x;
        particles[i].y = p_h[i].y + particles[0].y;
        particles[i].z = p_h[i].z + particles[0].z;
    }
}
