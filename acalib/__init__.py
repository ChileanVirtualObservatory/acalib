from . import core
from .core import *

<<<<<<< HEAD
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
=======
from . import algorithms

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore', category=UserWarning, append=True)
>>>>>>> f212972b038469446ffa48c5d21e86f18c311105

#from __future__ import absolute_import, print_function

#__all__ = []	

#from . import core
#from .core import *

#__all__.extend(core.__all__)
