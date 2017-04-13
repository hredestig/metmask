from __future__ import print_function
from builtins import str
from builtins import object
import xml.dom.minidom

import metmask.query as mmquery
from metmask.mask import mask
from metmask.parse.main import parserError


class chemspiderParser(object):
    def __init__(self, parent):
        self.token = parent.token
        if not self.token:
            raise parserError("""a security token is needed for using chemspider. 
register a user at http://www.chemspider.com and put your security token in the 
metmask configuration file in 'general' section under the name 'token'
""")
        search = 'http://www.chemspider.com/Search.asmx/'
        self.SimpleSearchURL = search + 'SimpleSearch?'
        self.GetCompoundInfoURL = search + 'GetCompoundInfo?'
        self.GetCompoundInfoURL = search + 'GetCompoundInfo?'
        self.ExtRefsURL = search + 'CSID2ExtRefs?'

        self.InChI = 'http://www.chemspider.com/InChI.asmx/'
        self.queryTables = ['smiles', 'chemspider', 'inchi', 'inchikey', 'synonym']
        self.parent = parent

    def InChItool(self, start, goal, query):
        urls = {'inchikey': {'chemspider': '/InChIKeyToCSID?inchi_key=',
                             'inchi': '/InChIKeyToInChI?inchi_key='},
                'inchi': {'chemspider': '/InChIToCSID?inchi=',
                          'inchikey': '/InChIToInChIKey?inchi=',
                          'smiles': '/InChIToSMILES?inchi='},
                'smiles': {'inchi': '/SMILESToInChI?smiles='}}
        if start not in urls:
            raise Exception("unsupported query")
        if goal not in urls[start]:
            raise Exception("unsupported query")
        query = self.parent.urlSafe(query)
        url = self.InChI + urls[start][goal] + query
        qRes = self.parent.getUrl(url)
        if not qRes:
            return ([])
        searchResults = xml.dom.minidom.parse(qRes)
        ids = list(mmquery.nodecontents(searchResults.getElementsByTagName('string')))
        return (ids)

    def ExtRefs(self, csid, ds):
        """
        take a csid
        """
        url = self.ExtRefsURL + "CSID=" + str(csid) + "&datasources=" + ds + "&token=" + self.token
        qRes = self.parent.getUrl(url)
        if not qRes:
            return ([])
        searchResults = xml.dom.minidom.parse(qRes)
        ids = list(mmquery.nodecontents(searchResults.getElementsByTagName('ext_id')))
        return (ids)

    def SimpleSearch(self, query):
        """
        take a query string, get the csids
        """
        query = self.parent.urlSafe(query)
        url = self.SimpleSearchURL + "query=" + str(query) + "&token=" + self.token
        qRes = self.parent.getUrl(url)
        if not qRes:
            return ([])
        searchResults = xml.dom.minidom.parse(qRes)
        ids = list(mmquery.nodecontents(searchResults.getElementsByTagName('int')))
        return (ids)

    def GetCompoundInfo(self, mm, csid):
        """
        take an csid, get a mask
        parameters:
        -`mm`: a metmask database
        -`csid`: a chemspider identifier
        """
        tmpmask = mask({}, mm.idpatterns)
        csid = self.parent.urlSafe(csid)
        url = self.GetCompoundInfoURL + "CSID=" + str(csid) + "&token=" + self.token
        qRes = self.parent.getUrl(url)
        if not qRes:
            return (tmpmask)
        tmpmask.append('chemspider', csid, self.parent.confid, self.parent.sourceid)
        searchResults = xml.dom.minidom.parse(qRes)
        ids = \
            list(mmquery.nodecontents(searchResults.getElementsByTagName('CompoundInfo')))
        inchi = list(mmquery.nodecontents(searchResults.getElementsByTagName('InChI')))
        if inchi:
            for ide in inchi:
                tmpmask.append('inchi', ide, \
                               self.parent.confid, self.parent.sourceid)
        smiles = \
            list(mmquery.nodecontents(searchResults.getElementsByTagName('SMILES')))
        if smiles:
            for ide in smiles:
                tmpmask.append('smiles', ide, self.parent.confid, self.parent.sourceid)
        inchikey = \
            list(mmquery.nodecontents(searchResults.getElementsByTagName('InChIKey')))
        if inchikey:
            for ide in inchikey:
                tmpmask.append('inchikey', ide, self.parent.confid, self.parent.sourceid)
        return (tmpmask)

    def getCsMasks(self, un, mm):
        """ take a mask and get filled mask from cs """
        res = {}
        tmpmask = mask({}, mm.idpatterns)
        if not any([x in self.queryTables for x in un.getTables()]):
            return (res)
        if un.hasTable('chemspider'):
            for csid in un.getIdentifiers('chemspider'):
                res[csid] = self.GetCompoundInfo(mm, csid)
                tmpmask.merge(res[csid])

        def filteredids(table, un, tmpmask):
            identifiers = []
            if un.hasTable(table):
                identifiers = un.getIdentifiers(table)
                if tmpmask.hasTable(table):
                    identifiers = [x for x in identifiers if x \
                                                   not in tmpmask.getIdentifiers(table)]
            return (identifiers)

        def fillresandmerge(ide, tmpmask, res):
            newids = self.SimpleSearch(ide)
            for csid in newids:
                if csid not in res:
                    res[csid] = self.GetCompoundInfo(mm, csid)
                    tmpmask.merge(res[csid])

        for ide in filteredids('smiles', un, tmpmask):
            fillresandmerge(ide, tmpmask, res)
        return (res)

    def fixInchi(self, un, assid=0):
        edited = False

        def addids(identifiers, table, assid):
            edited = False
            for ide in identifiers:
                if un.hasTable(table):
                    if ide in un.getIdentifiers(table):
                        continue
                un.append(table, ide, self.parent.confid, self.parent.sourceid, assid)
                edited = True
            return (edited)

        if assid == 0:
            assid = self.parent.mm.addAss()

        if un.hasTable('inchi'):
            for ide in un.getIdentifiers('inchi'):
                edited = addids(self.InChItool('inchi', 'chemspider', ide), \
                                'chemspider', assid)
                edited = addids(self.InChItool('inchi', 'inchikey', ide), \
                                'inchikey', assid)
                edited = addids(self.InChItool('inchi', 'smiles', ide), \
                                'smiles', assid)
        if un.hasTable('smiles'):
            for ide in un.getIdentifiers('smiles'):
                edited = addids(self.InChItool('smiles', 'inchi', ide), \
                                'inchi', assid)
        if un.hasTable('inchikey'):
            for ide in un.getIdentifiers('inchikey'):
                edited = addids(self.InChItool('inchikey', 'chemspider', ide), \
                                'inchikey', assid)
                edited = addids(self.InChItool('inchikey', 'inchi', ide), \
                                'inchikey', assid)
        if edited:
            self.fixInchi(un, assid)
        return (edited)


class parser(object):
    def __init__(self, parent):
        self.cs = chemspiderParser(parent)
        parent.tables = self.cs.queryTables
        parent.master = 'chemspider'
        self.parent = parent
        if not parent.boost:
            parent.boost = True

    def process(self):
        parent = self.parent
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
                if self.cs.fixInchi(un):
                    parent.setMask(un, setass=False)
                    tmp2 = parent.mm.getMask(tmp, weak=True)
                    if tmp2:
                        un = tmp2[0]
            try:
                csMasks = self.cs.getCsMasks(un, parent.mm)
            except:
                pdb.set_trace()
            if parent.mm.debug:
                print("#COMMENT mask " + str(ll) + " csid " + str(list(csMasks.keys())))
            cmsk = mask({})
            for csid in list(csMasks.keys()):
                csMasks[csid].setAllAssoc(parent.mm.addAss())
                cmsk.merge(csMasks[csid])
            parent.setMask(cmsk, setass=False)
