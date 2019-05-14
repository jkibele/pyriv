#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import os
from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

#this is the list of dependencies that pip uses to install
#requirements = [
#    'numpy',
#    'click',
#    'networkx',
#]

#extra requirements for optional features with their own dependencies
#on install, include the named extras in square brackets after the project name..
#i.e. pip install pyriv [geoFunc]
#Note: can also specify this with entry points for dynamic loading..
#i.e. ... ['pyriv=pyriv.cli:main [geoFunc]'], ['pyriv=pyriv.cli:other [otherExtra]'], ...
extra_requirements = {
    'geoFunc': ["geopandas"],
}

setup_requirements = [
    # TODO(jkibele): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pyriv',
    version='0.2.10',
    description="pyriv calculates minimum aquatic distance between points.",
    long_description=readme + '\n\n' + history,
    author="Jared Kibele",
    author_email='jkibele@nceas.ucsb.edu',
    url='https://github.com/jkibele/pyriv',
    download_url='https://github.com/jkibele/pyriv/releases',
    packages=find_packages(include=['pyriv'], exclude=['docs']),
    entry_points={
        'console_scripts': [
            'pyriv=pyriv.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=['numpy','networkx'],#requirements,
    extras_require=extra_requirements,
    license="MIT license",
    zip_safe=False,
    keywords='river distance pyriv',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: GIS',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
