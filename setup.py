# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('agatta/__version__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(

    name='agatta',

    version=main_ns['__version__'],

    python_requires='>=3.7.0',

    packages=find_packages(),

    author="Valentin Rineau",

    author_email="valentin.rineau@gmail.com",

    description="Three-item analysis python package",

    long_description=open('README.md').read(),

    install_requires= ["docopt",
                       "ete3",
                       "tqdm",
                       "treeswift"],

    include_package_data=True,

    url='https://github.com/vrineau/agatta',

    entry_points = {
        'console_scripts': [
            'agatta = agatta:main',
        ],
    },

    license="GNU GPLv3"

)

