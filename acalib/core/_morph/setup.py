from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
     ext_modules=[Extension("morph", ["morphWrapper.c", "morph.c"], extra_compile_args = ["-std=c11"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
