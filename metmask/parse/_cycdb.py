import pdb
import re

import metmask.parse
from metmask.mask import mask


class parser :

    def __init__ (self, parent) :
        self.uniqueid = parent.source + 'id'
        parent.tables = [self.uniqueid, 'cas', 'synonym', 'smiles', 'inchi', 'kegg']
        self.parent = parent

    def process (self) :
        """ parse the compunds.dat file from the Cyc databases and insert the
        cycdb frame-id (uniqueid), synonyms, kegg and cas identifiers
        """
        parent = self.parent
        un = mask({}, parent.mm.idpatterns) # new mask
        un.MIN_OVERLAP = 1
        ll = parent.getLine(comment='#')
        while ll:
            ll = ll.strip()
            if ll == '//' :
                if un.nid() > 2 :
                    try:
                        parent.setMask(un)
                    except:
                        #parent.mm.setMask(un, debug=True)
                        raise
                un = mask({}, parent.mm.idpatterns)
                un.MIN_OVERLAP = 1
            # frameid
            x = re.findall("UNIQUE-ID - (.*)", ll)
            if x :
                un.append(self.uniqueid, x[0], parent.confid, \
                              parent.sourceid)

            # synonym
            x = re.findall("COMMON-NAME - (.*)", ll) 
            if x :
                un.append('synonym', x[0], \
                              parent.mm.confidence['weak'], \
                              parent.sourceid)
            x = re.findall("SYNONYMS - (.*)", ll) 
            if x :
                un.append('synonym', x[0], \
                              parent.mm.confidence['weak'], \
                              parent.sourceid)

            # cas
            x = re.findall("DBLINKS - \(CAS \"(\d{1,7}-\d{2}-\d{1})\"", ll) 
            if x :
                un.append('cas', x[0], parent.confid, \
                              parent.sourceid)

            #inchi
            x = re.findall("INCHI - (.*)", ll) 
            if x :
               un.append('inchi', x[0], parent.confid, \
                             parent.sourceid)

            #smiles
            x = re.findall("SMILES - (.*)", ll) 
            if x :
               un.append('smiles', x[0], parent.mm.confidence['weak'], \
                             parent.sourceid)

            # kegg
            x = re.findall("DBLINKS - \(LIGAND-CPD \"([A-Z][0-9]{5})\"", ll) 
            if x :
                un.append('kegg', x[0], parent.confid, \
                              parent.sourceid)
            ll = parent.getLine(comment='#')
        if un.nid() > 2 :
            parent.setMask(un)
