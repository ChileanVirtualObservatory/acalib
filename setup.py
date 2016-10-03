from distutils.core import run_setup
from setuptools import setup, find_packages
from setuptools.command.install import install
from shutil import copyfile

import os
import glob

# PyCupid Library command
class BuildCupid(install):
    def run(self):
        cwd = os.getcwd()
        rel_pycupid = 'acalib/cupid/'
        pycupid_dir = os.path.join(cwd, rel_pycupid)
        os.chdir(pycupid_dir)

        run_setup('setup.py', ["build"])

        for fbuilded in glob.glob("build/lib*/*.so"):
            dest_directory = os.getcwd() + '/' + os.path.basename(fbuilded)
            copyfile(fbuilded, dest_directory)

        os.chdir(cwd)

setup(
    name = "acalib",

    version = "1.0.0",
    description = "Advanced Computing for Astronomy Library",
    url = "https://github.com/ChileanVirtualObservatory/ACALIB",
    author = "LIRAE",
    author_email = 'contact@lirae.cl',
    cmdclass = {"bcupid": BuildCupid},
    classifiers = [
        'Intended Audience :: Science/Research',

        'Topic :: Scientific/Engineering :: Astronomy',

        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'
        ],
    
    packages = find_packages(),
    include_package_data = True,
    setup_requires = ['numpy>=1.11', 'cython>=0.18'],
    install_requires = ['numpy>=1.11', 'astropy>=1.2', 'cython>=0.24',
                        'matplotlib>=1.5', 'scipy>=0.18',
                        'scikit-image>=0.12']
)
