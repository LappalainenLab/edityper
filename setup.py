#!/usr/bin/env python

"""Setup EdiTyper"""

from __future__ import division
from __future__ import print_function

#   Deal with Tkinter stuff
#   Needed for matplotlib
#   However, not on Pip
import sys
# try:
#     if sys.version_info.major == 2:
#         import Tkinter
#     elif sys.version_info.major == 3:
#         import tkinter
#     else:
#         sys.exit("Unsupported Python version")
# except ImportError as error:
#     sys.exit("Please install Tkinter and Tcl/Tk: " + str(error))


import os

#   Get stuff from setuptools
from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.install import install as _install

#   Load which
from edityper.toolkit import which

#   Some basic information
NAME = 'EdiTyper'
VERSION = '1.0.0'
AUTHOR = 'Alexandre Yahi'
AUTHOR_EMAIL = 'ay2318@cumc.columbia.edu'
DESCRIPTION = ''
LICENSE = ''
KEYWORDS = 'crispr rnaseq'
URL = 'https://github.com/lappalainenlab/EdiTyper'

#   Ensure Cython is available
if 'Cython.Distutils' not in sys.modules:
    try:
        import pip
    except ImportError:
        sys.exit("Please install Cython before installing %s" % NAME)
    INSTALL_CYTHON = ['install', 'cython']
    DEFAULT_DIR = tuple(filter(os.path.isdir, sys.path))[0]
    if not os.access(os.path.join(DEFAULT_DIR, 'site_packages'), os.W_OK):
        INSTALL_CYTHON.insert(1, '--user')
    pip.main(INSTALL_CYTHON)


from Cython.Distutils import build_ext

# File finder for R files
def find_r_files(dirname=None): # type: (Optional[str]) -> List[str]
    """Find R files"""
    rfiles = list()
    if not dirname:
        dirname = '.'
    for root, _, filenames in os.walk(dirname):
        for fname in filenames:
            if os.path.splitext(fname)[-1].upper() == '.R':
                rfiles.append(os.path.join(root, fname))
    return rfiles


#   A class to force install to run build first
class install(_install):
    """Force install to build first"""

    def run(self):
        self.run_command('build_ext')
        try:
            rscript = which('Rscript')
        except ValueError:
            raise
        else:
            pass
        _install.run(self)


#   Classifiers
CLASSIFIERS = [
    #   See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    #   For more classfiiers
    #   How mature is this project? Common values are
    'Development Status :: 5 - Production/Stable',
    #   What environment does this run in?
    'Environment :: Console',
    #   Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Visualization',
    #   Language
    'Natural Language :: English',
    #   Pick your license as you wish (should match "license" above)
    # 'License :: OSI Approved :: MIT License',
    #   Specify the Python versions you support here. In particular, ensure
    #   that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: C++',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: Implementation :: CPython',
    #   Operating systems we support
    # 'Operating System :: Microsoft :: Windows',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
]

#   Specify Python version
#   We support Python 2.7 and 3.5 or higher
PYTHON_REQUIRES='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4'

#   Platforms
# PLATFORMS = ['Windows', 'Linux', 'Mac OS-X', 'UNIX']
PLATFORMS = ['Linux', 'Mac OS-X', 'UNIX']

#   Dependencies
INSTALL_REQUIRES = [
    'cython',
    'numpy',
    'scipy',
    'matplotlib',
    'biopython',
    'regex'
]

#   Packages
PACKAGE_DIR = 'edityper'
PACKAGES = [
    PACKAGE_DIR
]

#   Set up the extension
EXT_MODULES = [
    Extension(
        name='%s.recnw' % PACKAGE_DIR,
        sources=['src/recnw_cython.pyx', 'src/recnw.cpp'],
        language='c++'
    )
]

#   Commands available for setup.py
CMD_CLASS = {
    'build_ext': build_ext,
    'install': install
}

#   Entry points into the program
ENTRY_POINTS = {
    'console_scripts': [
        '%(name)s = %(package)s.%(package)s:main' % {'name': NAME, 'package': PACKAGE_DIR}
    ]
}

#   Package data (R scripts)
PACKAGE_DATA = {'': ['*.r', '*.R']}

#   Run setup
setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    platforms=PLATFORMS,
    python_requires=PYTHON_REQUIRES,
    classifiers=CLASSIFIERS,
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,
    package_data=PACKAGE_DATA,
    include_package_data=True,
    ext_modules=EXT_MODULES,
    entry_points=ENTRY_POINTS,
    cmdclass=CMD_CLASS
)
