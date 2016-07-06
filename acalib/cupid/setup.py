from distutils.core import setup
from distutils.extension import Extension
import numpy as np

sourcefiles = []
wrapper_sources= ['mers.c','ast.c','pycf.c', 'cf.c', ]
cupidsub_sources = ['cupidcfaddpixel.c', 'cupidcfclump.c','cupidcfdeleteps.c',
'cupidcferode.c','cupidcfidl.c','cupidcfnebs.c','cupidcfscan.c',
'cupidcfxfer.c','cupidcfmakeps.c','cupidcffreeps.c','cupiddefminpix.c',
'cupidconfigrms.c','cupidcflevels.c','cupidconfigI.c','cupidconfigD.c']

sourcefiles += wrapper_sources
cupidsub_routes = ['cupidsub/'+c for c in cupidsub_sources]
sourcefiles += cupidsub_routes

setup(
  name = 'Cupid Library for Python',
  include_dirs = [np.get_include(),'include'],       
  ext_modules = [Extension("cupidpy",sourcefiles )]
)
