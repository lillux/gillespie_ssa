#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 11:44:53 2023

@author: lillo
"""
from setuptools import setup, find_packages

setup(name = 'gillespie_ssa',
     version = '0.1',
     description = 'Gillespie Stochastic Simulation Algorithm Library',
     url = 'https://github.com/lillux/gillespie_ssa',
     author = 'Calogero Carlino',
     author_email = 'calogero.carlino28@gmail.com',
     license = 'GPLv3',
     zip_safe=False,
     install_requires=['numpy>=1.18.1',
                       'matplotlib'],
                       
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.8',
          'Topic :: Software Development :: Libraries',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
              
      packages=find_packages())