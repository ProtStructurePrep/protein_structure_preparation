# -*- coding: utf-8 -*-
"""
Setup for repo template.

The project contains a **single top-level package**.
The top-level package package.py stores all the project info.
"""

import codecs
from os import path

from setuptools import find_packages, setup


def get_long_description(readme):
    """Get the long description from the README file."""
    with codecs.open(readme, encoding='utf-8') as _rf:
        return _rf.read()


def get_package_name():
    'The top-level package name.'
    top_level_packages = [
        p for p in find_packages(exclude=['tests']) if '.' not in p
    ]
    if len(top_level_packages) != 1:
        raise ValueError('Project must contain a single top-level package.')
    return top_level_packages[0]


base_dir = path.abspath(path.dirname(__file__))
readme_file = path.join(base_dir, 'README.md')

# single top-level package
package_name = get_package_name()

# get project info from the package.py module of the top-level package
project_info = {}
with open(path.join(package_name, 'package.py')) as fp:
    exec(fp.read(), project_info)  # pylint: disable=exec-used

project_name = project_info['__title__']
version = project_info['__version__']

long_description = get_long_description(readme_file)

packages = find_packages(exclude=['tests'])
modules = []

SETUP_REQUIRES = []
INSTALL_REQUIRES = [
    # # this is an example of URL based requirement (see PEP508):
    # 'repo @ http://github.com/user/repo/archive/master.tar.gz',
]
EXTRAS_REQUIRES = {}

entry_points = project_info.get('__entry_points__')
if entry_points is not None:
    entry_points = [f'{key}={val}' for key, val in entry_points.items()]

setup(
    name=project_name,
    version=version,
    description=project_info.get('__summary__'),
    long_description=long_description,
    author=project_info.get('__author__'),
    author_email=project_info.get('__email__'),
    url=project_info.get('__url__'),
    py_modules=modules,
    packages=packages,
    setup_requires=SETUP_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRES,
    license=project_info.get('__license__'),
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=project_info.get('__classifiers__'),
    entry_points={
        'console_scripts': entry_points,
    },
    scripts=project_info.get('__scripts__'))
