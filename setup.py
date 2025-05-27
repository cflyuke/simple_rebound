from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import numpy

# Define the C extension modules
ext_modules = [
    Extension(
        'simple_rebound.core',
        sources=[
            os.path.join('src', 'core_module.c'),
            os.path.join('src', 'collision.c'),
            os.path.join('src', 'gravity.c'),
            os.path.join('src', 'integrator.c'),
            os.path.join('src', 'integrator_ias15.c'),
            os.path.join('src', 'integrator_mercurius.c'),
            os.path.join('src', 'integrator_whfast.c'),
            os.path.join('src', 'particle.c'),
            os.path.join('src', 'simulation.c'),
            os.path.join('src', 'transformations.c'),
        ],
        include_dirs=['src', numpy.get_include()],
        extra_compile_args=['-O3'],
    )
]

# Custom build command to handle platform-specific compilation
class CustomBuildExt(build_ext):
    def build_extensions(self):
        if os.name == 'nt':  # Windows
            for ext in self.extensions:
                ext.extra_compile_args = ['-O2']
        super().build_extensions()

setup(
    name="simple_rebound",
    version="0.1.0",
    author="cflyuke",
    description="A simplified N-body simulation package",
    packages=["simple_rebound"],
    ext_modules=ext_modules,
    cmdclass={'build_ext': CustomBuildExt},
    python_requires=">=3.6",
    zip_safe=False,
)
