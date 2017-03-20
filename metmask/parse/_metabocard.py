from metmask.mask import mask
import re,pdb
import metmask.parse
from main import parserError
from metmask.mask import idMisMatchError

class parser:

    def __init__(self, parent) :
        parent.tables = ['biocyc', 'cas', 'chebi', 'formula', 'hmdb', 'inchi', \
                             'iupac','kegg', 'metlin','synonym', 'sid', 'cid', \
                             'smiles',]
        self.parent = parent
        self.tableDict = {'# biocyc_id:': [parent.confid, 'biocyc'],
                          '# cas_number:': [parent.confid, 'cas'],
                          '# chebi_id:' : [parent.confid, 'chebi'],
                          '# chemical_formula:': [parent.mm.confidence['weak'], 'formula'],
                          '# hmdb_id:': [parent.confid, 'hmdb'],
                          '# inchi_identifier:': [parent.confid, 'inchi'],
                          '# iupac:': [parent.confid, 'iupac'],
                          '# kegg_compound_id:': [parent.confid, 'kegg'],
                          '# metlin_id:': [parent.confid, 'metlin'],
                          '# name:' : [parent.confid, 'synonym'],
                          '# synonyms:' : [parent.confid, 'synonym'],
                          '# pubchem_compound_id:' : [parent.confid, 'cid'],
                          '# pubchem_substance_id:' : [parent.confid, 'sid'], 
                          '# smiles_isomeric:' : [parent.confid, 'smiles']}

    def process (self) :
        """ parse a file with metabocards from hmdb 
        """
        parent = self.parent
        un = mask({}, parent.mm.idpatterns) 
        ll = parent.getLine(comment='\n')
        if not ll.startswith("#BEGIN_METABOCARD") :
            raise parserError, "file does seem to contain metabocards: " + ll
        parent.mm.setTableWeak('formula')
        ok = False

        while ll :
            if ll.startswith('#END_METABOCARD') :
                if un.nid(parent.mm) > 1 :
                    parent.setMask(un)

            elif ll.startswith("#BEGIN_METABOCARD") :
                un = mask({}, parent.mm.idpatterns)

            elif ll.startswith("#") :
                ok = False
                if self.tableDict.has_key(ll.strip()) :
                    con = self.tableDict[ll.strip()][0]
                    tab = self.tableDict[ll.strip()][1]
                    ok = True

            else :
                if ok :
                    identifiers = ll.split("; ")
                    for iden in identifiers :
                        if not iden.lower().startswith("not available"):
                            try :
                                un.append(tab, iden, con, parent.sourceid)
                            except idMisMatchError:
                                print "#OFFENDING LINE " + \
                                    str(parent.lineNum) + "@" + \
                                    tab + " : " + str(iden)
            ll = parent.getLine(comment='\n')

      
