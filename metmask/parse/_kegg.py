import re

from .main import fileFormatError
from metmask.mask import mask


class parser:
    """ Parse and import the PATHWAYS from the KEGG compounds database
    dump file. KEGG identifiers are created for those compounds that
    aren't already in the database but no other identifiers are read.
    """

    def __init__(self, parent):
        self.parent = parent
        parent.tables = ['kegg', 'synonym', 'cas', 'formula',
                         'chebi', 'knapsack', 'sid', 'pathway']
        parent.master = 'kegg'

    def process(self):
        parent = self.parent
        # annotate this table as fully weak
        parent.mm.setTableWeak('pathway')
        parent.mm.setTableWeak('formula')

        ll = parent.getLine()
        if not re.match("^ENTRY", ll):
            raise fileFormatError("COMPOUNDS file doesn't look like expected")

        un = mask({}, parent.mm.idpatterns)
        good = True

        while ll:
            # check and send off gathered data
            # kegg
            x = re.findall('ENTRY *([A-Z][0-9]{5}) ', ll)
            if x:
                if not un.isEmpty() and good:
                    parent.setMask(un)

                un = mask({}, parent.mm.idpatterns)  # new mask
                keggno = x[0].strip()
                un.append('kegg', keggno, parent.confid,
                          parent.sourceid)
                good = True

            # synonym
            x = re.findall('NAME +(.+)', ll)
            if x:
                un.append('synonym',
                          x[0].strip().replace(";", ""),
                          parent.confid, parent.sourceid)
                while re.match(".+;$", ll):
                    ll = parent.getLine()
                    un.append('synonym',
                              ll.strip().replace(";", ""),
                              parent.confid,
                              parent.sourceid)

            # formula
            x = re.findall('FORMULA +(.+)', ll)
            if x:
                un.append('formula', x[0], parent.mm.confidence['weak'],
                          parent.sourceid)
            # knapsack
            x = re.findall('KNApSAcK: (.+)', ll)
            if x:
                un.append('knapsack', x[0], parent.confid,
                          parent.sourceid)
            # cas
            if re.findall('CAS:', ll):
                x = re.findall('(\d{2,7}-\d{2}-\d{1})', ll)
                for xx in x:
                    un.append('cas', xx, parent.confid,
                              parent.sourceid)
            # sid <- obs not cid!
            x = re.findall('PubChem: (\d+)', ll)
            if x:
                un.append('sid', x[0], parent.confid,
                          parent.sourceid)
            # chebi
            x = re.findall('ChEBI: (\d+)', ll)
            if x:
                un.append('chebi', x[0], parent.confid,
                          parent.sourceid)
            # pathway
            x = re.findall('(PATHWAY)*\s*(ko|map)(\d+)', ll)
            if x:
                un.append('pathway', x[0][2],
                          parent.mm.confidence['weak'],
                          parent.sourceid)

            # comment, to get rid of kegg internal generics
            if re.match('ENTRY *[A-Z][0-9]{5}$ *Peptide *Compound$', ll):
                good = False
            if re.match('COMMENT', ll):
                if re.match('COMMENT     generic compound in reaction hierarchy', ll) or \
                        re.match('COMMENT     coordination compound', ll) or \
                        re.match('COMMENT.+R=.+', ll):
                    good = False

            ll = parent.getLine()

        # send of the last mask if we have one
        if not un.isEmpty() and good:
            parent.setMask(un)
