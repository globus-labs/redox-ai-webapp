from setuptools import setup, find_packages

# this grabs the requirements from requirements.txt
requirements = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name='redoxweb',
    version='0.0.1',
    packages=find_packages(),
    description='Web service for designing molecules online',
    install_requires=requirements
)
