from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
import numpy as np


cythonize('morph.pyx')

setup(
	name= 'Morphological Operators for Spectra Sketcher',
    ext_modules = [Extension("morph",sources=["morph.c", "c_morph.c"])],
	include_dirs=[np.get_include()],
)
