from setuptools import setup, find_packages

setup(
    name='igh_converter',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
    ],
    author='marietai',
    description='A package for converting SNP positions to a custom IGH reference.',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
