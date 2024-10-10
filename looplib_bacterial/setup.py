import numpy as np
from setuptools import setup
from Cython.Build import cythonize

# The setup script
setup(
    name="looplib_bacterial",
    version="0.1",
    install_requires=['numpy'],
    description="looplib by Anton Goloborodko, edited for circular bacterial chromosomes",
    packages=["looplib_bacterial"],
    ext_modules=cythonize(['looplib_bacterial/bacterial_no_bypassing.pyx','looplib_bacterial/bacterial_bypassing.pyx']),
    include_dirs=[np.get_include()]
)
