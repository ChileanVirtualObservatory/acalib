import sys
sys.path.append('../../')

from acalib.core import *
from os import walk
from astropy import log

binpath='../../bindata/fits/'


f = []
for (dirpath, dirnames, filenames) in walk(binpath):
   for filename in filenames:
      f.append(dirpath+'/'+filename)
cs =[]
for fn in filter(lambda x: '.fits' in x, f):
   c = AContainer()
   log.info("Loading "+fn)
   c.load_from_fits(fn)
   cs.append(c)





