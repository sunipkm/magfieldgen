from setuptools import setup #, Extension
import os, os.path
import re

longDescription= "3D Stochastic Magnetic Field Generator using a C backend, which in turn uses FFTW3 and OpenMP multithreading for speed up. Dynamic Creator Mersenne Twister is used for thread safe random number generation. Hence, FFTW3, and DCMT are necessary for successful build and execution of this package."


setup(name='magfieldgen',
      version='2.0',
      description='3D Stochastic Magnetic Field Generator using C libraries',
      author='Sunip Mukherjee',
      author_email='sunipkmukherjee@gmail.com',
      license='New BSD',
      long_description=longDescription,
      package_dir = {'magfieldgen/': ''},
      packages=['magfieldgen'],
      install_requires=['numpy']
      )
