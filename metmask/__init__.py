"""
A package for managing metabolite identifiers for Metabolomics
experiments.

Copyright (C) Henning Redestig
2009
See COPYING.txt for licensing details.
"""

__all__ = ['dbi', 'parse', 'mask', 'query']
__version__ = '0.5.4'

WEAK_CONF = 2
""" reserved confidence index for weak associations  """
NEVERMERGE_CONF = 1
""" reserved confidence index for nevermerge associations  """
