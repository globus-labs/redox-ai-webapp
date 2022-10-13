from setuptools import setup, find_packages

setup(
    name='redoxweb',
    version='0.0.1',
    packages=find_packages(),
    description='Web service for designing molecules online',
    install_requires=[
        'rdkit-pypi',
        'fastapi',
        'dlhub-sdk',
        'foundry_ml'
    ]
)
