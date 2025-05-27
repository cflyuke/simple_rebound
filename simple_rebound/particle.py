"""
Particle class for Simple Rebound
"""

import numpy as np

class Particle:
    """
    表示一个粒子的Python类
    """
    
    def __init__(self, m=0.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, r=0.0):
        """
        初始化粒子
        
        参数:
            m: 质量
            x, y, z: 位置
            vx, vy, vz: 速度
            r: 半径
        """
        self.m = float(m)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.vx = float(vx)
        self.vy = float(vy)
        self.vz = float(vz)
        self.r = float(r)
        self.ax = 0.0
        self.ay = 0.0
        self.az = 0.0
        self.id = -1
    
    @property
    def position(self):
        """返回位置向量"""
        return np.array([self.x, self.y, self.z])
    
    @position.setter
    def position(self, pos):
        """设置位置向量"""
        self.x, self.y, self.z = float(pos[0]), float(pos[1]), float(pos[2])
    
    @property
    def velocity(self):
        """返回速度向量"""
        return np.array([self.vx, self.vy, self.vz])
    
    @velocity.setter
    def velocity(self, vel):
        """设置速度向量"""
        self.vx, self.vy, self.vz = float(vel[0]), float(vel[1]), float(vel[2])
    
    @property
    def acceleration(self):
        """返回加速度向量"""
        return np.array([self.ax, self.ay, self.az])
    
    @property
    def speed(self):
        """返回速度大小"""
        return np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)
    
    @property
    def kinetic_energy(self):
        """返回动能"""
        return 0.5 * self.m * (self.vx**2 + self.vy**2 + self.vz**2)
    
    @property
    def momentum(self):
        """返回动量向量"""
        return self.m * np.array([self.vx, self.vy, self.vz])
    
    def distance_to(self, other):
        """计算到另一个粒子的距离"""
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return np.sqrt(dx**2 + dy**2 + dz**2)
    
    def relative_velocity_to(self, other):
        """计算相对于另一个粒子的相对速度"""
        dvx = self.vx - other.vx
        dvy = self.vy - other.vy
        dvz = self.vz - other.vz
        return np.sqrt(dvx**2 + dvy**2 + dvz**2)
    
    def set_orbit(self, central_mass, a, e=0.0, inc=0.0, Omega=0.0, omega=0.0, M=0.0, G=6.67430e-11):
        """
        设置轨道参数
        
        参数:
            central_mass: 中心天体质量
            a: 半长轴
            e: 偏心率
            inc: 倾角
            Omega: 升交点经度
            omega: 近心点幅角
            M: 平近点角
            G: 引力常数
        """
        # 简化实现：假设圆轨道
        r = a
        v = np.sqrt(G * central_mass / r)
        
        # 设置位置和速度
        self.x = r * np.cos(M)
        self.y = r * np.sin(M)
        self.z = 0.0
        
        self.vx = -v * np.sin(M)
        self.vy = v * np.cos(M)
        self.vz = 0.0
        
        # 应用倾角旋转
        if inc != 0.0:
            cos_i = np.cos(inc)
            sin_i = np.sin(inc)
            
            new_y = self.y * cos_i - self.z * sin_i
            new_z = self.y * sin_i + self.z * cos_i
            self.y = new_y
            self.z = new_z
            
            new_vy = self.vy * cos_i - self.vz * sin_i
            new_vz = self.vy * sin_i + self.vz * cos_i
            self.vy = new_vy
            self.vz = new_vz
    
    def copy(self):
        """创建粒子的副本"""
        p = Particle(self.m, self.x, self.y, self.z, self.vx, self.vy, self.vz, self.r)
        p.ax = self.ax
        p.ay = self.ay
        p.az = self.az
        p.id = self.id
        return p
    
    def __str__(self):
        """字符串表示"""
        return (f"Particle(m={self.m:.3e}, "
                f"pos=({self.x:.3e}, {self.y:.3e}, {self.z:.3e}), "
                f"vel=({self.vx:.3e}, {self.vy:.3e}, {self.vz:.3e}), "
                f"r={self.r:.3e})")
    
    def __repr__(self):
        """详细字符串表示"""
        return self.__str__()
