from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import re

import metmask
import metmask.parse
from .main import fileFormatError, fixLine
from metmask.mask import idMisMatchError, mask


class parser(object):
    """ Parse and import identifiers from a file which is formatted
    so
    that the header names the table name (exactly, e.g. to have
    something imported as a CAS number it must be under the header
    ``cas`` and nothing else.) and that the string underneath, *the
    major field*, is a list of identifiers. The list is one one or
    more identifiers separated by the specified charactor, *the minor
    field*. Separators characters are allowed inside the major field
    but not the minor. Returns the number of successfully processed
    masks.

    A table name can be preceeded by 'weak' to make that table hold
    annotations rather than identifiers.

    Each line is interpreted as a single mask.
    All lines must have the
    same number of major fields. It is preferable to quote all
    fields.

    Example::

      "cas","synonym"\n
      "id1|id2","n_id1|n_id2"\n
    """

    def __init__(self, parent):
        tabs = fixLine(parent.getLine(), sep=parent.sep1)
        self.confids = list(range(0, len(tabs)))
        for i in range(0, len(tabs)):
            x = re.findall('weak:(.+)', tabs[i])
            self.confids[i] = parent.confid
            if x:
                tabs[i] = x[0].strip()
                self.confids[i] = metmask.WEAK_CONF
                # parent.mm.setTableWeak(tabs[i])
        parent.tables = tabs
        self.parent = parent

    def process(self):
        parent = self.parent
        ncol = len(parent.tables)
        ## single identifier inserts should be allowed since they can
        ## might map to bins
        # if ncol < 2:
        #    raise fileFormatError, "Only one column, pointless insertion" 
        # make sure that the necessary tables are in the db
        ll = parent.getLine()
        ll = ll.replace("#", "\\#")
        while ll:
            ll = ll.replace("#", "\\#")
            try:
                vec = fixLine(ll, sep=parent.sep1)
                if len(vec) != ncol:
                    raise fileFormatError("Number of columns doesn't match the header :" \
                        + str(parent.lineNum) + ll)
                un = mask({}, parent.mm.idpatterns)
                for i in range(0, len(vec)):
                    idvec = vec[i].strip().split(parent.sep2)
                    idvec[0].strip()
                    if idvec[0] and not re.match(parent.na, idvec[0]):
                        for ide in idvec:
                            ide = ide.strip()
                            if not re.match(parent.na, ide):
                                try:
                                    un.append(parent.tables[i], ide,
                                              self.confids[i],
                                              parent.sourceid)
                                except idMisMatchError:
                                    print("#OFFENDING LINE " + \
                                          str(parent.lineNum) + "@" + \
                                          parent.tables[i] + \
                                          " : " + str(ide))
                # no empty masks
                if not un.isEmpty():
                    parent.setMask(un)
            except KeyboardInterrupt:
                raise Exception("Interrupt caught, breaking")
            except fileFormatError:
                print("#ERROR: format problem" + ll)
            except:
                print("#ERROR:" + ll)
            ll = parent.getLine()
