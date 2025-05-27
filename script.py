#!/usr/bin/env python3
"""
测试简化版rebound包的示例
重现原始代码的功能
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random
import math
import time

# 导入我们的简化rebound包
import simple_rebound.simple_rebound as simple_rebound

# 设置随机种子
seed = 42
random.seed(seed)
np.random.seed(seed)

class Config():
    def __init__(self):
        """
        实验参数设定：
        1. 自变量配置
        2. 实验配置：实验时长控制，记录频率
        """
        ## 一、自变量配置
        # G值
        self.G = 6.67430e-11
        # 边界（正方体边长）
        self.boundary_size = 8e11
        # 碰撞处理机制：‘merge’代表完全非弹性碰撞
        self.collision = 'merge'
        # 碰撞检测机制：direct(default), spatial
        self.collision_detection = 'direct'
        # 积分器类型：ias15(default), leapfrog, whfast, mercurius
        self.integrator = "mercurius"
        # 重力计算机制：direct(default)，tree
        self.gravity_method = "direct"
        # 积分步长(s)
        self.dt = 86400.0

        self.R_mean = 1.5e11  # 轨道半径(1Au)
        self.M_star = 1.989e30 # 恒星质量
        self.r = 6.96e8 #恒星半径
        self.N_dust = 100 # 总星子数
        self.rho = 3000  # 密度 kg/m³
        self.M_total = 6e24  # 星子总质量（这里是地球）

        self.R_std = self.R_mean / 200  # 初始半径扰动项
        self.Z_std = 0 # 初始z轴坐标标准差

        self.V_mean = math.sqrt(self.G * self.M_star / self.R_mean) 
        self.V_std = self.V_mean / 200 # 初始x,y合速度标准差
        self.VZ_std = 0 # 初始z轴速度标准差
        
        self.M_mean = self.M_total / self.N_dust  # 质量采用gauss分布，改变分布需要同时改变下面的质量设定代码
        self.M_std = self.M_mean / 10 

        ## 二、实验配置
        self.update_step = 86400 * 365 * 100  # 信息更新记录时间(s)
        self.total_time = 86400 * 365 * 1000   # 总体时间(s)


        


def create_terrestrial_system(config):
    """创建类地行星形成系统"""

    
    # 创建模拟
    sim = simple_rebound.Simulation()
    sim.G = config.G
    sim.boundary_size = config.boundary_size
    sim.collision = config.collision
    sim.collision_detection = config.collision_detection
    sim.integrator = config.integrator
    sim.gravity_method = config.gravity_method
    sim.dt = config.dt
    R_mean = config.R_mean
    M_star = config.M_star
    r = config.r
    sim.add(m=M_star, r=r)
    N_dust = config.N_dust
    rho = config.rho
    M_total = config.M_total
    R_std = config.R_std
    Z_std = config.Z_std
    V_mean = config.V_mean
    V_std = config.V_std
    VZ_std = config.VZ_std 
    M_mean = M_total / N_dust
    M_std = config.M_std
    
    def rotate(theta, x, y):
        """旋转坐标"""
        return (x * math.cos(theta) - y * math.sin(theta), 
                x * math.sin(theta) + y * math.cos(theta))
    
    # 添加尘埃粒子
    for i in range(N_dust):
        r = random.gauss(R_mean, R_std) 
        x1, y1, z = r, 0, random.gauss(0, Z_std)
        vx1, vy1, vz = (random.gauss(0, V_std), 
                        random.gauss(V_mean, V_std), 
                        random.gauss(0, VZ_std))
        
        theta = random.random() * 2 * math.pi
        x, y = rotate(theta, x1, y1)
        vx, vy = rotate(theta, vx1, vy1)
        
        m = random.gauss(M_mean, M_std) # 质量设定
        radius = (m / rho * 3 / (4 * math.pi)) ** (1/3)
        
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, r=radius)
    
    return sim

def realtime_visualization():
    """实时可视化模拟过程"""
    print("Simple Rebound real time simulation begin:")
    print("=" * 50)
    
    # 创建模拟
    config = Config()
    sim= create_terrestrial_system(config)
    print(f"Original particle number: {sim.N-1}")
    plt.style.use('seaborn-v0_8-darkgrid')
    colors = plt.cm.tab20.colors
    # 设置交互模式
    plt.ion()
    fig = plt.figure(figsize=(12, 12))
    
    # 3D轨迹图
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-2, 2)
    ax1.set_zlim(-0.5, 0.5)
    ax1.set_xlabel('X [AU]')
    ax1.set_ylabel('Y [AU]')
    ax1.set_zlabel('Z [AU]')
    ax1.set_title('Particle Trajectories')
    ax1.view_init(elev=30, azim=45)
    
    # 初始散点图
    particles_scatter = ax1.scatter([], [], [], s=10, c='blue', alpha=0.6, depthshade=True)
    time_text = ax1.text2D(0.02, 0.95, '', transform=ax1.transAxes)
    
    # 模拟参数
    update_step = config.update_step
    total_time = config.total_time
    target_fps = 10  # 目标帧率
    min_frame_time = 1.0 / target_fps  # 最小帧间隔(秒)
    
    step_count = 0
    last_frame_time = time.time()
    max_mass = []
    min_mass = []
    particle_number = []
    while sim.t < total_time:
        # 积分一步
        sim.integrate(sim.t + update_step)
        step_count += 1
        
        # 更新可视化前控制帧率
        current_time = time.time()
        elapsed = current_time - last_frame_time
        if elapsed < min_frame_time:
            time.sleep(min_frame_time - elapsed)
        last_frame_time = time.time()
        
        # 记录初始质量分布(第一次迭代时)
        if step_count == 1:
            initial_masses = [p.m for p in sim.particles[1:]]
        
        # 获取当前粒子状态(忽略太阳)
        particles = []
        for p in sim.particles[1:]:
            particles.append({
                'x': p.x, 'y': p.y, 'z': p.z,
                'vx': p.vx, 'vy': p.vy, 'vz': p.vz,
                'm': p.m, 'r': p.r
            })
        particle_number.append(sim.N-1)
        particle_m = np.array([p['m'] for p in particles])
        max_mass.append(np.max(particle_m))
        min_mass.append(np.min(particle_m))
        
        if particles:
            # 更新3D散点图
            x = [p['x'] / 1.5e11 for p in particles]
            y = [p['y'] / 1.5e11 for p in particles]
            z = [p['z'] / 1.5e11 for p in particles]
            particles_scatter._offsets3d = (x, y, z)
            
            # 更新大小和颜色(根据质量着色)
            sizes = np.array([math.pow(p['m'] / config.M_mean, 1/3) * 5 for p in particles], dtype=np.float32)
            sizes = np.maximum(sizes, 1)
            particles_scatter.set_sizes(sizes)
            
            # 使用质量对数作为颜色值
            masses = np.array([p['m'] for p in particles], dtype=np.float32)
            log_masses = np.log10(masses)
            particles_scatter.set_array(log_masses)
            particles_scatter.set_cmap('viridis')  # 使用viridis色图
            
            # 更新时间显示
            current_time = sim.t / (86400 * 365)
            time_text.set_text(f'Time: {current_time:.1f} years\nParticles: {sim.N-1}')
            
            # 重绘并确保帧率
            fig.canvas.draw()
            fig.canvas.flush_events()
        
        # 如果粒子数太少就停止
        if sim.N <= 2:
            print(f"The number of particle is {sim.N -1}, simulation shuts down.")
            break
    
    print(f"Simulation complete, final particle number: {sim.N-1}")
    plt.ioff()

    # 保存动画为GIF格式
    ani = FuncAnimation(fig, lambda i: None, frames=len(max_mass), interval=100)
    ani.save('simulation.gif', writer='pillow', fps=10)
    # 保存质量变化曲线图
    fig2, ax = plt.subplots(figsize=(10, 6))
    ax.plot(max_mass, color=colors[0], label='Max Mass')
    ax.plot(min_mass, color=colors[1], label='Min Mass')
    ax.set_title('Mass Evolution Over Time')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Mass (kg)')
    ax.legend()
    plt.savefig('mass_evolution.png')
    plt.close(fig2)

    # 保存粒子数目变化图
    fig3, ax = plt.subplots(figsize=(10, 6))
    ax.plot(particle_number, color=colors[2], label='Particle Number')
    ax.set_title('Particle Number Evolution')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Number of Particles')
    ax.legend()
    plt.savefig('particle_number_evolution.png')
    plt.close(fig3)
    
    # 绘制初始和最终质量分布对比
    plt.figure(figsize=(12, 6))
    plt.subplot(121)
    plt.hist(initial_masses, bins=20, alpha=0.7, color=colors[3], label='Initial')
    plt.title('Initial Mass Distribution')
    plt.xlabel('Mass (kg)')
    plt.ylabel('Count')
    
    plt.subplot(122)
    final_masses = [p['m'] for p in particles]
    plt.hist(final_masses, bins=20, alpha=0.7, color=colors[4], label='Final')
    plt.title('Final Mass Distribution')
    plt.xlabel('Mass (kg)')
    plt.ylabel('Count')
    
    plt.tight_layout()
    plt.savefig('mass_distribution_comparison.png')
    plt.show()

def main():
    """主函数"""
    realtime_visualization()
    print("Completed!")

if __name__ == "__main__":
    main()
