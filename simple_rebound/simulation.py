"""
Simulation class for Simple Rebound - C backend only
"""

import numpy as np
from .particle import Particle
from . import core

class Simulation:
    """
    N体模拟的主要类 - 使用C扩展后端
    """
    
    def __init__(self):
        """初始化模拟"""
        self._c_sim = core.Simulation()
    
    def add(self, m=0.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, r=0.0, particle=None):
        """
        添加粒子到模拟中
        
        参数:
            m: 质量
            x, y, z: 位置
            vx, vy, vz: 速度
            r: 半径
            particle: Particle对象（如果提供，其他参数将被忽略）
        """
        if particle is not None:
            m, x, y, z = particle.m, particle.x, particle.y, particle.z
            vx, vy, vz, r = particle.vx, particle.vy, particle.vz, particle.r
        
        self._c_sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, r=r)
    
    def remove(self, index):
        """移除指定索引的粒子"""
        self._c_sim.remove(index)
    
    @property
    def particles(self):
        """获取粒子列表（从C扩展动态生成）"""
        c_particles = self._c_sim.particles
        python_particles = []
        
        for i, c_p in enumerate(c_particles):
            p = Particle(c_p['m'], c_p['x'], c_p['y'], c_p['z'],
                        c_p['vx'], c_p['vy'], c_p['vz'], c_p['r'])
            p.id = c_p['id']
            python_particles.append(p)
        
        return python_particles
    
    def get_particle(self, index):
        """获取指定索引的粒子"""
        c_p = self._c_sim.get_particle(index)
        p = Particle(c_p['m'], c_p['x'], c_p['y'], c_p['z'],
                    c_p['vx'], c_p['vy'], c_p['vz'], c_p['r'])
        p.id = c_p['id']
        return p
    
    @property
    def N(self):
        """返回粒子数量"""
        return self._c_sim.N
    
    @property
    def t(self):
        """获取当前时间"""
        return self._c_sim.t
    
    @t.setter
    def t(self, value):
        """设置当前时间"""
        self._c_sim.t = value
    
    @property
    def dt(self):
        """获取时间步长"""
        return self._c_sim.dt
    
    @dt.setter
    def dt(self, value):
        """设置时间步长"""
        self._c_sim.dt = value
    
    @property
    def G(self):
        """获取引力常数"""
        return self._c_sim.G
    
    @G.setter
    def G(self, value):
        """设置引力常数"""
        self._c_sim.G = value
    
    @property
    def softening(self):
        """获取软化系数"""
        return self._c_sim.softening
    
    @softening.setter
    def softening(self, value):
        """设置软化系数"""
        self._c_sim.softening = value
    
    @property
    def theta(self):
        """获取Barnes-Hut开角参数"""
        return self._c_sim.theta
    
    @theta.setter
    def theta(self, value):
        """设置Barnes-Hut开角参数"""
        self._c_sim.theta = value
    
    @property
    def gravity_method(self):
        """获取引力计算方法"""
        return self._c_sim.gravity_method
    
    @gravity_method.setter
    def gravity_method(self, value):
        """设置引力计算方法"""
        self._c_sim.gravity_method = value
    
    @property
    def integrator(self):
        """获取积分器类型"""
        return self._c_sim.integrator
    
    @integrator.setter
    def integrator(self, value):
        """设置积分器类型"""
        self._c_sim.integrator = value
    
    @property
    def collision(self):
        """获取碰撞类型"""
        return self._c_sim.collision
    
    @collision.setter
    def collision(self, value):
        """设置碰撞类型"""
        self._c_sim.collision = value
    
    @property
    def collision_detection(self):
        """获取碰撞检测方法"""
        return self._c_sim.collision_detection
    
    @collision_detection.setter
    def collision_detection(self, value):
        """设置碰撞检测方法"""
        self._c_sim.collision_detection = value
    
    @property
    def boundary(self):
        """获取边界类型"""
        return self._c_sim.boundary
    
    @boundary.setter
    def boundary(self, value):
        """设置边界类型"""
        self._c_sim.boundary = value
    
    @property
    def boundary_size(self):
        """获取边界大小"""
        return self._c_sim.boundary_size
    
    @boundary_size.setter
    def boundary_size(self, value):
        """设置边界大小"""
        self._c_sim.boundary_size = float(value)
    
    def set_boundary(self, boundary="open", size=100.0):
        """
        设置边界条件
        
        参数:
            boundary: 边界类型 ("open", "periodic", "reflective")
            size: 边界大小（立方体边长）
        """
        self._c_sim.set_boundary(boundary=boundary, size=size)
    
    def step(self):
        """执行一个时间步"""
        self._c_sim.step()
    
    def integrate(self, tmax, exact_finish_time=True):
        """
        积分到指定时间
        
        参数:
            tmax: 目标时间
            exact_finish_time: 是否精确到达目标时间
        """
        if exact_finish_time:
            # 精确到达目标时间的积分
            while self.t < tmax:
                # 检查下一步是否会超过目标时间
                if self.t + self.dt > tmax:
                    # 计算到达目标时间所需的精确步长
                    remaining_time = tmax - self.t
                    old_dt = self.dt
                    self.dt = remaining_time
                    # 对于IAS15等自适应积分器，需要临时禁用自适应步长
                    self._c_sim.set_exact_finish_time(True)
                    self.step()
                    self.dt = old_dt
                    self._c_sim.set_exact_finish_time(False)
                    break
                else:
                    self.step()
        else:
            # 不需要精确到达，使用标准积分
            while self.t < tmax:
                self.step()
    
    def calculate_energy(self):
        """计算系统总能量"""
        return self._c_sim.calculate_energy()
    
    def calculate_momentum(self):
        """计算系统总动量"""
        return self._c_sim.calculate_momentum()
    
    def status(self):
        """打印模拟状态"""
        print(f"Simulation Status (C extension):")
        print(f"  Time: {self.t:.6e}")
        print(f"  Particles: {self.N}")
        print(f"  Time step: {self.dt:.6e}")
        print(f"  G: {self.G:.6e}")
        print(f"  Softening: {self.softening:.6e}")
        print(f"  Theta: {self.theta:.6e}")
        print(f"  Integrator: {self.integrator}")
        print(f"  Collision: {self.collision}")
        print(f"  Collision detection: {self.collision_detection}")
        print(f"  Gravity method: {self.gravity_method}")
        
        # 边界信息
        boundary_str = self.boundary
        if self.boundary != "open":
            boundary_str += f" (size: {self.boundary_size:.6e})"
        print(f"  Boundary: {boundary_str}")
        
        print(f"  Energy: {self.calculate_energy():.6e}")
        momentum = self.calculate_momentum()
        print(f"  Momentum: ({momentum[0]:.6e}, {momentum[1]:.6e}, {momentum[2]:.6e})")
