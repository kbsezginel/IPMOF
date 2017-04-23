import os
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="ipmof",
    version="0.1",
    description="Interpenetrating MOFs",
    author="Kutay B. Sezginel",
    author_email="kbs37@pitt.edu",
    install_requires=requirements,
    packages=find_packages(),
    include_package_data=True
)
