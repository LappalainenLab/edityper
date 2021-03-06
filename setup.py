#!/usr/bin/env python

"""Setup EdiTyper"""

from __future__ import division
from __future__ import print_function

import os
import sys

#   Builtins to tell EdiTyper if we're setting up or not
#   Inspired by NumPy
if sys.version_info.major == 3:
    import builtins
else:
    import __builtin__ as builtins


builtins.__EDITYPER_SETUP__ = True

#   Get stuff from setuptools
from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.install import install as _install

#   Some basic information
NAME = 'EdiTyper'
AUTHOR = 'Alexandre Yahi'
AUTHOR_EMAIL = 'ay2318@cumc.columbia.edu'
DESCRIPTION = ''
LICENSE = 'MIT License'
LICENSE_FILE = 'LICENSE'
KEYWORDS = 'crispr rnaseq'
URL = 'https://github.com/LappalainenLab/edityper'

#   Ensure Cython is available
try:
    import pip
except ImportError:
    sys.exit("Please install Cython before installing %s" % NAME)

if 'cython' not in {mod.key for mod in pip.get_installed_distributions()}:
    INSTALL_CYTHON = ['install', 'cython']
    DEFAULT_DIR = tuple(filter(os.path.isdir, sys.path))[0]
    if not os.access(os.path.join(DEFAULT_DIR, 'site_packages'), os.W_OK):
        INSTALL_CYTHON.insert(1, '--user')
    pip.main(INSTALL_CYTHON)


from Cython.Distutils import build_ext

from edityper import __version__ as VERSION

#   A class to force install to run build first
class install(_install):
    """Force install to build first"""

    def run(self):
        self.run_command('build_ext')
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
    'License :: OSI Approved :: MIT License',
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
    license_file=LICENSE_FILE,
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
