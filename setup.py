#!/usr/bin/env python

from setuptools.core import setup
import os
import metmask
import metmask.parse

files = ["data/*"]

setup(name='metmask',
      version=metmask.__version__,
      description='A program for masking metabolite identifiers',
      author='Henning Redestig',
      author_email='henning.red@gmail.com',
      url='http://metmask.sourceforge.net',
      requires=['sqlite3', 'SOAPpy'],
      platforms=['Linux', 'WinXP'],
      classifiers = [
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bioinformatics',
        ],
      license='OSI Approved :: GNU General Public License (GPL)',
      packages=['metmask', 'metmask.parse'],
      long_description="""
      This is a package for creating, maintaining and querying a
      database with metabolite identifiers. Focused on mapping analyte
      identifiers to the original identifiers of the parent metabolite
      in order to facilitate biological interpretation of metabolomics
      datasets. Provides automated import for several sources such as
      KEGG, PlantCyc and the NIST library.
      """,
      package_data={'metmask': files},
      scripts=['scripts/metmask'])


