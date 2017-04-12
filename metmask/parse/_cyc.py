import re

from .main import fixLine, parserError
from metmask.mask import mask


class parser:
    def __init__(self, parent):
        parent.tables = ['synonym', 'cas', 'kegg', 'cycpath', 'smiles',
                         'formula', 'pubchem', 'knapsack']
        parent.sep1 = "\t"
        header = parent.getLine()
        vec = fixLine(header, sep=parent.sep1)
        self.parent = parent
        try:
            vec = [vec[0]] + [vec[1]] + [vec[5]]
            if vec != ['Compound_common_name',
                       'Compound_synonyms', 'Links']:
                raise parserError("file does not look like expected")
        except:
            raise parserError("file does not look like expected:" + header)

    def process(self):
        """ parse a cyc dump file and fectch cas, kegg, synonym and
        cycpath
        """
        parent = self.parent
        last = ""
        parent.mm.setTableWeak('cycpath')
        parent.mm.setTableWeak('formula')

        l = parent.getLine()
        while l:
            l = l.replace("'", "\\'")  # cyc likes apostrophes
            l = l.replace("#", "\#")  # hashes in the middle of lines
            vec = fixLine(l, sep=parent.sep1)

            if len(vec) < 6:
                l = parent.getLine()
                continue
            if len(vec) == 9:
                vec = [vec[0]] + [vec[1]] + [vec[5]] + [vec[8].strip()] + [vec[3]] + [vec[4]]
            else:
                vec = [vec[0]] + [vec[1]] + [vec[5]] + [''] + [vec[3]] + [vec[4]]

            if vec == last:
                l = parent.getLine()
                continue
            un = mask({}, parent.mm.idpatterns)
            if len(vec) == 6:
                # pointless if we dont have links
                if vec[2] != '':
                    # cas
                    x = re.findall('CAS:(\d{1,7}-\d{2}-\d{1})',
                                   vec[2])
                    for casno in x:
                        un.append('cas', casno, parent.confid,
                                  parent.sourceid)
                    # kegg
                    x = re.findall('LIGAND-CPD:([A-Z][0-9]{5})',
                                   vec[2])
                    for kegg in x:
                        un.append('kegg', kegg, parent.confid,
                                  parent.sourceid)
                    # pubchem
                    x = re.findall('PUBCHEM:([0-9]+)',
                                   vec[2])
                    for kegg in x:
                        un.append('pubchem', kegg, parent.confid,
                                  parent.sourceid)
                    # knapsack
                    x = re.findall('KNAPSACK:([0-9]+)',
                                   vec[2])
                    for kegg in x:
                        un.append('knapsack', kegg, parent.confid,
                                  parent.sourceid)
                    # formula
                    if vec[4] != '':
                        un.append('formula', vec[4].replace(" ", ""),
                                  parent.mm.confidence['weak'],
                                  parent.sourceid)
                    # smiles
                    if vec[5] != '':
                        un.append('smiles', vec[5],
                                  parent.mm.confidence['weak'],
                                  parent.sourceid)
                    # synonym
                    un.append('synonym', vec[0], parent.mm.confidence['weak'],
                              parent.sourceid)
                    if vec[1] != '':
                        x = vec[1].split("*")
                        for syn in x:
                            un.append('synonym', syn,
                                      parent.mm.confidence['weak'],
                                      parent.sourceid)
                    if vec[3] != '':
                        un.append('cycpath', vec[3],
                                  parent.mm.confidence['weak'],
                                  parent.sourceid)
                    if un.nid() < 2 or \
                            all([not un.hasTable('cas'), not un.hasTable('kegg')]):
                        un = mask({}, parent.mm.idpatterns)
            if un.nid() > 1:
                parent.setMask(un)
            last = vec
            l = parent.getLine()
