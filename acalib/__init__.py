from . import core
from .core import *

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

#from __future__ import absolute_import, print_function

#__all__ = []	

#from . import core
#from .core import *

#__all__.extend(core.__all__)
