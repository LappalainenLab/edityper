#!/usr/bin/env python2

"""Setup CRISPronto"""

from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.install import install as _install
from Cython.Distutils import build_ext

#   Some basic information
NAME = 'CRISPronto'
VERSION = '1.0'
AUTHOR = 'Alexandre Yahi'
DESCRIPTION = ''
LICENSE = ''
KEYWORDS = ''
URL = ''

#   A class to force install to run build first
class install(_install):
    """Force install to build first"""

    def run(self):
        self.run_command('build_ext')
        _install.run(self)


#   Set up the extension
EXT_MODULES = [
    Extension(
        name='NW_py',
        sources=['NW_py.pyx', 'py_align.cpp'],
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
        '%s = test:main' % NAME
    ]
}

#   Run setup
setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    ext_modules=EXT_MODULES,
    entry_points=ENTRY_POINTS,
    cmdclass=CMD_CLASS
)
