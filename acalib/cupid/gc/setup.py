from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

sources = ["pycupid.pyx"]

wrapper_sources = [
#"mers.c",
#"msgOutiff.c",
#"msgOutif.c", 
#"msgBlankif.c", 
#"emsMark.c",
#"ems1Mpush.c",
#"ems1Emark.c"
]

cupid_sources = [
#"cupidgaussclumps.c", 
#"cupidconfigrms.c", 
#"cupiddefminpix.c", 
#"cupidedges.c",
#"cupidgccalcf.c", 
#"cupidgccalcg.c",
#"cupidgcchisq.c",
#"cupidgcdump.c", 
#"cupidgcfindmax.c",
#"cupidgcfit.c",
#"cupidgclistclump.c",
#"cupidgcmodel.c",
#"cupidgcndfclump.c",
#"cupidgcprofwidth.c", 
#"cupidgcsetinit.c",
#"cupidgcupdatearrays.c",
#"cupidndfclump.c",  
#"cupidrca2.c",
#"cupidrca.c",
#"cupidrcopyline.c",
#"cupidconfigI.c",
#"cupidconfigD.c"
]

cupid_sources = ["cupidsub/" + s for s in cupid_sources]

sources += cupid_sources + wrapper_sources

extensions = [
    Extension(
        "pycupid",
        sources,
        include_dirs = ["./libraries/ast/include", "./include"],
        library_dirs = ["./libraries/ast/lib"],
        libraries = ["ast"],
        extra_link_args = [
        "-last",
        "-last_pal", 
        "-last_grf_2.0", 
        "-last_grf_3.2",
        "-last_grf_5.6",
        "-last_grf3d",
        "-last_pass2", 
        "-last_err", 
        "-lm",
        ],
    ),
]

setup(
    name = "PyCupid",
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(extensions),
)