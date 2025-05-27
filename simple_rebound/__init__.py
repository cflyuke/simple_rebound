"""
Simple Rebound - A simplified N-body simulation package
"""

from .simulation import Simulation
from .particle import Particle

__version__ = "0.1.0"
__author__ = "Simple Rebound Team"

# 导出主要类
__all__ = ['Simulation', 'Particle']
