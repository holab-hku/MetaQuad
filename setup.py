"""
MetaQuad, a pipeline for calling informative mutations in shotgun metagenomic data
See: https://github.com/SX-16/MetaQuad
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

local_path = path.abspath(path.dirname(__file__))
#local_path = os.path.join(path.abspath(path.dirname(__file__)), "metaquad")

# Set __version__
exec(open("metaquad/version.py").read())

#Description
with open(path.join(local_path, 'README.md'), encoding='utf-8') as file:
    read_me = file.read()

requirement = ['numpy>=1.9.0', 'matplotlib', 'BBMix>=0.2.1', 'vireoSNP', 'pandas', 'seaborn']

setup(
    name='metaquad',
    version=__version__,
    description='MetaQuad, a pipeline for calling informative mutations in shotgun metagenomic data',
    long_description=read_me,
    url='https://github.com/SX-16/MetaQuad',
    author='Sheng XU',
    author_email='shengxu@connect.hku.hk',
    license='',
    keywords=['informative mutation', 'shotgun metagenomic sequencing',
              'beta binomial mixture model'],


    packages=find_packages(),

    entry_points={
        'console_scripts': [
            'metaquad = metaquad.main:main',
        ],
    },


    install_requires=requirement,

    extras_require={
        'docs': [
            'sphinx_bootstrap_theme'
        ]
    },

    py_modules=['metaquad']

)
