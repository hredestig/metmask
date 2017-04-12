import re
import socket

import xml.dom.minidom

# import SOAPpy
# from SOAPpy import WSDL

import metmask.query as mmquery
from metmask.mask import guessTable, mask

socket.setdefaulttimeout(10)


class pubchemParser:
    def __init__(self, parent):
        self.wsdlserver = \
            WSDL.Proxy('http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/eutils.wsdl')
        self.queryTables = ['iupac', 'cid', 'smiles', 'inchi', 'inchikey', 'formula',
                            'weight', 'totalcharge', 'xlogp', 'hbonddonor',
                            'hbondacceptor', 'heavyatom', 'tpsa']
        self.parent = parent
        self.MAX_RET = str(1)
        self.mail = parent.mm.MAIL
        self.tool = parent.mm.METMASK
        self.fields = {'cid': 'UID',
                       'inchi': 'INCHI',
                       'inchikey': 'INCHIKEY',
                       'iupac': 'IUPAC',
                       'kegg': 'SYNO',
                       'cas': 'SYNO'}
        self.db = 'pccompound'
        self.queryResult = {'WebEnv': '', 'QueryKey': ''}
        self.summary = \
            'http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=' + self.db + \
            '&tool=' + self.tool + \
            '&mail=' + self.mail + \
            '&retmax' + self.MAX_RET

    # def queryPubchem(self, un, table):
    #     """ perform queries for table. results are left online """
    #     if table not in self.fields:
    #         raise Exception("unsupported query")
    #     if un.hasTable(table):
    #         for ide in un.getIdentifiers(table):
    #             qres = self.wsdlserver.run_eSearch(tool=self.tool, email=self.mail, term="\"" + ide + "\"",
    #                                                field=self.fields[table], retmax=self.MAX_RET, db=self.db,
    #                                                usehistory='y')
    #             if isinstance(qres, SOAPpy.Types.structType):
    #                 if qres['Count'] != '0':
    #                     self.queryResult = qres
    #                     return (True)
    #     return (False)

    def pubchem2mask(self, docSum):
        """ turn a docsum node in to a mask
        """
        un = mask({}, self.parent.mm.idpatterns)
        cid = next(mmquery.nodecontents(docSum.getElementsByTagName("Id")))
        un.append('cid', cid, self.parent.confid, self.parent.sourceid)
        p = re.compile('(<a href[^>]*>)|(</a>)|(ligand)|(,)')
        cnf = self.parent.confid
        src = self.parent.sourceid
        weak = self.parent.mm.confidence['weak']

        for item in docSum.getElementsByTagName('Item'):
            if len(item.childNodes) == 0:
                continue
            nVal = item.childNodes[0].nodeValue
            if item.getAttribute('Name') == 'SynonymList':
                synonyms = list(mmquery.nodecontents(item.childNodes))
                for s in synonyms:
                    s = p.sub('', s.strip())
                    if 'kegg' in guessTable(s, self.parent.mm.idpatterns):
                        un.append('kegg', s, cnf, src)
                        continue
                    if 'cas' in guessTable(s, self.parent.mm.idpatterns):
                        un.append('cas', s, cnf, src)
                        continue
                    if 'chebi' in guessTable(s, self.parent.mm.idpatterns):
                        un.append('chebi', s, cnf, src)
                        continue
                    if 'inchi' in guessTable(s, self.parent.mm.idpatterns):
                        un.append('inchi', s, cnf, src)
                        continue
                    else:
                        un.append('synonym', s, weak, src)
            elif item.getAttribute('Name') == 'IUPACName':
                un.append('iupac', nVal, cnf, src)
            elif item.getAttribute('Name') == 'CanonicalSmile':
                un.append('smiles', nVal, cnf, src)
            elif item.getAttribute('Name') == 'CanonicalSmile':
                un.append('smiles', nVal, cnf, src)
            elif item.getAttribute('Name') == 'InChIKey':
                un.append('inchikey', nVal, cnf, src)

            # annotation section
            elif item.getAttribute('Name') == 'MolecularFormula':
                un.append('formula', nVal, weak, src)
            elif item.getAttribute('Name') == 'MolecularWeight':
                un.append('weight', nVal, weak, src)
            elif item.getAttribute('Name') == 'TotalFormalCharge':
                un.append('totalcharge', nVal, weak, src)
            elif item.getAttribute('Name') == 'XLogP':
                un.append('xlogp', nVal, weak, src)
            elif item.getAttribute('Name') == 'XLogP':
                un.append('xlogp', nVal, weak, src)
            elif item.getAttribute('Name') == 'HydrogenBondDonorCount':
                un.append('hbonddonor', nVal, weak, src)
            elif item.getAttribute('Name') == 'HydrogenBondAcceptorCount':
                un.append('hbondacceptor', nVal, weak, src)
            elif item.getAttribute('Name') == 'HeavyAtomCount':
                un.append('heavyatom', nVal, weak, src)
            elif item.getAttribute('Name') == 'TPSA':
                un.append('tpsa', nVal, weak, src)
        return (un)

    def getPubchemMasks(self, un):
        """ perform queries cid, inchi, inchikey in un. make masks of all results """
        hit = False
        res = []
        if self.queryPubchem(un, 'cid'):
            hit = True
        elif self.queryPubchem(un, 'inchi'):
            hit = True
        elif self.queryPubchem(un, 'inchikey'):
            hit = True
        elif self.queryPubchem(un, 'iupac'):
            hit = True
        elif self.queryPubchem(un, 'kegg'):
            hit = True
        elif self.queryPubchem(un, 'cas'):
            hit = True
        if not hit:
            return (res)
        try:
            summaryUrl = self.summary + '&WebEnv=' + self.queryResult['WebEnv'] + "&query_key=" + self.queryResult[
                'QueryKey']
            summaryResults = xml.dom.minidom.parse(self.parent.getUrl(summaryUrl)).documentElement
            candidates = summaryResults.getElementsByTagName('DocSum')
            for cand in candidates:
                pcmask = self.pubchem2mask(cand)
                res.append(pcmask)
                break
            return (res)
        except:
            return ([])


class parser:
    def __init__(self, parent):
        self.pc = pubchemParser(parent)
        parent.tables = self.pc.queryTables
        parent.master = 'cid'
        self.parent = parent
        if not parent.boost:
            parent.boost = True

    def process(self):
        parent = self.parent
        weakTables = ['formula', 'weight', 'totalcharge', 'xlogp', 'hbonddonor',
                      'hbondacceptor', 'heavyatom', 'tpsa']
        list(map(lambda x: parent.mm.setTableWeak(x), weakTables))
        ll = True

        while ll:
            ll = parent.getLine()
            if not ll:
                break
            tmp = mask({})
            tmp.append('_id', ll)
            tmp2 = parent.mm.getMask(tmp, weak=True)
            if not tmp2:
                continue
            else:
                un = tmp2[0]
            try:
                pcMasks = self.pc.getPubchemMasks(un)
            except:
                print("#COMMENT unknown error")
                continue
            if parent.mm.debug:
                print("#COMMENT mask " + str(ll))
            msk = mask({})
            for pcmsk in pcMasks:
                pcmsk.setAllAssoc(parent.mm.addAss())
                msk.merge(pcmsk)
            parent.setMask(msk, setass=False)
