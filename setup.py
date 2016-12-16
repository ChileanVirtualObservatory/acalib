from setuptools import setup, find_packages
from setuptools.command.install import install
from shutil import copyfile

import os
import glob

import sys
import subprocess 

def check_build():
    good_commands = ('develop', 'sdist', 'build', 'build_ext', 'build_py',
                     'build_clib', 'build_scripts', 'bdist_wheel', 'bdist_rpm',
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg')

    for command in good_commands:
        if command in sys.argv[1::]:
            return True

def build_and_move(path):
    cwd = os.getcwd()
    rel_module = path
    module_dir = os.path.join(cwd, rel_module)
    os.chdir(module_dir)
    subprocess.call(["python","setup.py", "build"])

    for fbuilded in glob.glob("build/lib*/*.so"):
        dest_directory = os.getcwd() + '/' + os.path.basename(fbuilded)
        copyfile(fbuilded, dest_directory)

    os.chdir(cwd)


def setup_package():
    #Building packages
    if check_build():
        build_and_move('acalib/cupid')
        build_and_move('acalib/core/_morph')

    setup(
        name = "acalib",
        version = "0.1.0",
        description = "Advanced Computing for Astronomy Library",
        url = "https://github.com/ChileanVirtualObservatory/ACALIB",
        author = "CSRG",
        author_email = 'contact@lirae.cl',
        classifiers = [
            'Intended Audience :: Science/Research',

            'Topic :: Scientific/Engineering :: Astronomy',

            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6'
            ],
        
        zip_safe = False,
        packages = find_packages(),
        include_package_data = True,
        setup_requires = ['numpy>=1.8', 'cython>=0.18'],
        install_requires = ['numpy>=1.8', 'astropy>=1.2', 'cython>=0.24',
                            'matplotlib>=1.5', 'scipy>=0.18',
                            'scikit-image>=0.12']
    )




setup_package()
