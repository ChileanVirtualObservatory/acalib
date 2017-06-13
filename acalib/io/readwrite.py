from astropy import log
from astropy.table import Table
from tqdm import tqdm
import tarfile
from acalib.io.fits import save_fits_from_cont,load_fits_to_cont, loadFITS_PrimaryOnly
import os
import urllib.request
import astropy.units as u
# TODO: Dummy function. This should be smart (Zip, Tar, etc)

def uncompress(file,dir=None):
    if dir is None:
        dir = ""
    tar = tarfile.open(file)
    tar.extractall(path=dir)
    tsize = []
    tname = []
    for f in tar.members:
        tname.append(f.name)
        tsize.append(f.size)
    T = Table()
    T['filename'] = tname
    T['size'] = (tsize * u.B).to("MB")
    tar.close()
    return T

# TODO: Dummy function. This should be smart (Tables, AData etc...)
def load(uri):
    return loadFITS_PrimaryOnly(uri)

# TODO: tqdm_notebook (under _tqdm_notebook.tqdm_notebook) does not show well. Using tqdm
class TqdmUpTo(tqdm):
    """Alternative Class-based version of the above.

    Provides `update_to(n)` which uses `tqdm.update(delta_n)`.

    Inspired by [twine#242](https://github.com/pypa/twine/pull/242),
    [here](https://github.com/pypa/twine/commit/42e55e06).
    """
    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize

def download(url,dir=None,name=None):
    if name is None:
        name = url.rsplit('/', 1)[-1]
    if dir is None:
        dir = ""
    filename = os.path.join(dir, name)
    # Download the file if not found
    if not os.path.isfile(filename):
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=name) as t:  # all optional kwargs
            urllib.request.urlretrieve(url, filename=filename, reporthook=t.update_to, data=None)
    else:
        log.info("File "+name+" is already downloaded, skipping...")
    return filename



class Container:
    """
    Data structure that contains a list of AData and astropy tables.
    
    It can load from and save to FITS file format. To write a FITS, please put in "primary"
    the main object, and in nndata and table all the extensions.
    """
    def __init__(self): 
        self.primary = None
        """Primary object""" 
        self.images = [] 
        """List of NDData object""" 
        self.tables = []
        """List of astropy tables""" 

    def load_fits(self,path): 
        load_fits_to_cont(path,self)
    def save_fits(self,path): 
        save_fits_from_cont(path,self) 


def load_fits(path):
    """
    Load a FITS into a container.

    Parameters
    ----------
    path : str
        Path to FITS file in local disk.

    Returns
    -------
    result: :class:`~acalib.Container` with the FITS loaded.
    """
    cont=Container()
    cont.load_fits(path)
    return cont

def save_fits(cont,path):
    """
    Create a new FITS file from a container.

    Parameters
    ----------
    cont : :class:`~acalib.Container`

    path : str
    	Path to new FITS file to be created from the container.
    """
    cont.save_fits(path)
