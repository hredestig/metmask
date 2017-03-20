from metmask.mask import mask
import re,pdb
import metmask.parse

class parser:

    def __init__(self, parent) :
        parent.tables = ['nist', 'cas', 'synonym', 'formula']
        self.parent = parent

    def process (self) :
        """ parse an SDF file and import any NISTID, synonyms, CASNO
        and formula
        """
        parent = self.parent
        un = mask({}, parent.mm.idpatterns) 
        ll = parent.getLine()
        parent.mm.setTableWeak('formula')

        while ll :
            if ll.startswith('$$$$') :
                if un.nid(parent.mm) > 1 :
                    parent.setMask(un)
                un = mask({}, parent.mm.idpatterns)
            # synonym
            if ll.startswith('>  <NAME>') or ll.startswith('>  <SYNONYMS>'):
                tmp = parent.getLine().strip()
                while tmp != "" :
                    un.append('synonym', tmp, \
                                  parent.confid, parent.sourceid)
                    tmp = parent.getLine().strip()
            # nist
            if ll.startswith('>  <NISTNO>') :
                tmp = parent.getLine().strip()
                while tmp != "" :
                    un.append('nist', tmp, \
                                  parent.confid, parent.sourceid)
                    tmp = parent.getLine().strip()
            # cas
            if ll.startswith('>  <CASNO>') :
                tmp = parent.getLine().strip()
                while tmp != "" :
                    un.append('cas', tmp, \
                                  parent.confid, parent.sourceid)
                    tmp = parent.getLine().strip()
            # formula
            if ll.startswith('>  <FORMULA>') :
                tmp = parent.getLine().strip()
                while tmp != "" :
                    un.append('formula', tmp, \
                                  parent.mm.confidence['weak'], parent.sourceid)
                    tmp = parent.getLine().strip()
            ll = parent.getLine()
        if un.nid(parent.mm) > 1 :
            parent.setMask(un)
