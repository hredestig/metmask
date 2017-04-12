import os
import re
import metmask.parse.main

tmp = [x for x in os.listdir(__path__[0]) if re.match("_[^_]", x)] + ['main']

__all__ = list(set([re.sub("\.py.*$", "", x) for x in tmp]))

tmp = [x for x in __all__ if re.match("_[^_]", x)]
PARSERS = tmp
""" the current list of available parsers  """

from metmask.parse import *

del tmp
