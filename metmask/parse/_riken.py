import re

from metmask.mask import mask


class parser:
    def __init__(self, parent):
        parent.tables = ['riken', 'cas', 'synonym', 'kegg', 'formula', 'smiles']
        self.parent = parent

    def process(self):
        """ parse the RIKEN standards library file
        """
        parent = self.parent
        un = mask({}, parent.mm.idpatterns)
        ll = parent.getLine()
        parent.mm.setTableWeak('formula')

        while ll:
            x = re.findall('^Name: (.+)', ll)
            if x:
                if un.nid() > 1:
                    parent.setMask(un)
                un = mask({}, parent.mm.idpatterns)
                un.append('riken', x[0].strip(),
                          parent.confid, parent.sourceid)
            x = re.findall('CASNO: (\d{1,7}-\d{2}-\d{1})', ll)
            if x:
                un.append('cas', x[0].strip(),
                          parent.confid, parent.sourceid)
            x = re.findall('^KEGG: ([A-Z][0-9]{5})', ll)
            if x:
                un.append('kegg', x[0].strip(),
                          parent.confid, parent.sourceid)
            x = re.findall('^Formula: (.+)', ll)
            if x:
                un.append('formula', x[0].strip(),
                          parent.mm.confidence['weak'], parent.sourceid)
            x = re.findall('^SMILES: (.+)', ll)
            if x:
                un.append('smiles', x[0].strip(),
                          parent.mm.confidence['weak'], parent.sourceid)
            x = re.findall('Synonym: Name:(.+)', ll)
            if x:
                un.append('synonym', x[0].strip(),
                          parent.confid, parent.sourceid)
            ll = parent.getLine()
        if un.nid() > 1:
            parent.setMask(un)
