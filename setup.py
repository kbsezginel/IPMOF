import os
from setuptools import setup, find_packages


setup(
    name="ipmof",
    version="0.1",
    description="Interpenetrating MOFs",
    author="Kutay B. Sezginel",
    author_email="kbs37@pitt.edu",
    install_requires=[
        'xlrd',
        'pyyaml',
        'numpy',
        'matplotlib',
        'ase',
        'tabulate'
    ],
    packages=find_packages(),
    include_package_data=True
)
