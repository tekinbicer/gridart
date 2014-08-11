#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension                

ext_gridart = Extension(name='gridart.libgridart',
                        sources=['gridart/art.c'])

# Main setup configuration.
setup(
      name='gridart',
      version=1,
      ext_modules=[ext_gridart],
      author='Doga Gursoy',
      author_email='dgursoy@aps.anl.gov',
      description='Grid computing for ART',
      keywords=['tomography', 'reconstruction', 'imaging'],
      license='BSD',
      platforms='Any',
      classifiers=['Development Status :: 4 - Beta',
		   'License :: OSI Approved :: BSD License',
		   'Intended Audience :: Science/Research',
		   'Intended Audience :: Education',
		   'Intended Audience :: Developers',
		   'Natural Language :: English',
		   'Operating System :: OS Independent',
		   'Programming Language :: Python',
		   'Programming Language :: Python :: 2.6',
		   'Programming Language :: Python :: 2.7',
		   'Programming Language :: C',
		   'Programming Language :: C++']
      )



