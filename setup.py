# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('agatta/__version__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(

    name='Agatta',

    version=main_ns['__version__'],

    python_requires='>=3.8.0',

    packages=find_packages(),

    author="Valentin Rineau",

    author_email="valentin.rineau@sorbonne-universite.fr",

    description="Three-item analysis python package",

    long_description_content_type="text/markdown",

    long_description=open('README.md').read(),

    install_requires= ["docopt==0.6.2",
                       "ete3==3.1.3",
                       "tqdm==4.66.3",
                       "treeswift==1.1.37",
                       "six==1.16.0",
                       "numpy==1.24.4",
                       "pandas==1.4.4",
                       "pypdf==3.17.0"],

    include_package_data=True,

    url='https://github.com/vrineau/Agatta',

    entry_points = {
        'console_scripts': [
            'agatta = agatta:main',
        ],
    },

    license="GNU General Public License v3.0+",

    classifiers = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: MacOS",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Education",
    "Intended Audience :: Other Audience",
    "Natural Language :: English",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
]

)

