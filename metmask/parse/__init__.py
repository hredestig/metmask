import os
import re
import metmask.parse.main

tmp = filter(lambda x: re.match("_[^_]", x), os.listdir(__path__[0])) + ['main']

__all__ = list(set(map(lambda x: re.sub("\.py.*$", "", x), tmp)))

tmp = filter(lambda x: re.match("_[^_]", x), __all__)
PARSERS = tmp
""" the current list of available parsers  """

from metmask.parse import *

del tmp
