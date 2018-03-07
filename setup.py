#!/usr/bin/env python

from setuptools import setup

setup(
    name='hiplexpipe_somatic',
    version='0.1',
    author='Khalid Mahmood',
    author_email='khalid.mahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['hiplexpipe_somatic = src.main:main']
    },
    url='https://github.com/khalidm/hiplexpipe_somatic',
    license='LICENSE.txt',
    description='hiplexpipe_somatic is a bioinformatics pipeline to call somatic variants from tumour normal HiPlex data.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
