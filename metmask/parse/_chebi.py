import re, xml.dom.minidom, urllib2, pdb, httplib

from metmask.mask import idMisMatchError
from metmask.mask import mask
import re,pdb
import metmask.parse 
import metmask.query as query
import socket

socket.setdefaulttimeout(5)

class chebiParser:
    
    def __init__ (self, parent) :
        self.getComplete = 'http://www.ebi.ac.uk/webservices/chebi/test/getCompleteEntity?'
        self.getLite = 'http://www.ebi.ac.uk/webservices/chebi/test/getLiteEntity?'
        self.getOntologyChildren = 'http://www.ebi.ac.uk/webservices/chebi/test/getOntologyChildren?'
        self.queryTables = ['chebi', 'iupac', 'cas', 'kegg', 'inchi',\
                                'smiles', 'formula', 'synonym']
        self.parent = parent

    def url2ids (self, url):
        """
        take an url, get the chebi ids
        """
        qRes = self.parent.getUrl(url)
        if qRes :
            searchResults = xml.dom.minidom.parse(qRes)
            ids = list(query.nodecontents(searchResults.getElementsByTagName('ns1:chebiId')))
            return(map(lambda x:x.replace('CHEBI:', ''), ids))
        else:
            return([])

    def getChebiChildren (self, chebiId) :
        """ get the is_enantiomer children 
        """
        url = self.getOntologyChildren + "chebiId=" + chebiId
        qRes = self.parent.getUrl(url)
        if not qRes :
            return([])
        searchResults = xml.dom.minidom.parse(qRes)
        res = []
        for it in searchResults.getElementsByTagName('ns1:ListElement') :
            type = list(query.nodecontents(it.getElementsByTagName("ns1:type")))
            if type[0] == 'is enantiomer of' :
                child = list(query.nodecontents(it.getElementsByTagName("ns1:chebiId")))
                res.append(child[0].replace("CHEBI:", ""))
        return(res)
        
    def getChebiMasks (self, un, mm):
        """
        1. use getLite to query for all relevant entries in the mask
        2. Query for chebi, iupac, cas, kegg, inchi, smiles
           1. Get only exact matches
        3. Take the unique ones
        4. Get all is_a children of the found chebis and add them
        Result is a list of masks
        """
        res = {}
        #pdb.set_trace()
        tmpmask = mask({}, mm.idpatterns)
        if not any(map(lambda x: x in self.queryTables, un.getTables())) :
            return(res)
        if un.hasTable('chebi') :
            for ch in un.getIdentifiers('chebi') :
                res[ch] = self.chebi2mask(mm, ch)
                tmpmask.merge(res[ch])

        def filteredids (table, un, tmpmask):
            identifiers = []
            if un.hasTable(table) :
                identifiers = un.getIdentifiers(table)
                if tmpmask.hasTable(table) :
                    # to start querying for the same identifiers several times
                    identifiers = filter(lambda x:x not in tmpmask.getIdentifiers(table),\
                                             identifiers)
            return(identifiers)

        def fillresandmerge (qUrl, tmpmask, res) :
            newids = self.url2ids(qUrl)
            for ch in newids :
                if not res.has_key(ch):
                    res[ch] = self.chebi2mask(mm, ch)
                    tmpmask.merge(res[ch])

        for ide in filteredids('iupac', un, tmpmask) :
            qUrl = self.getLite + "search=" + \
                self.parent.urlSafe(ide) + "&searchCategory=IUPAC+NAME"
            fillresandmerge(qUrl, tmpmask, res)

        for ide in filteredids('cas', un, tmpmask) :
            qUrl = self.getLite + "search=" + \
                self.parent.urlSafe(ide) + "&searchCategory=REGISTRY+NUMBER"
            fillresandmerge(qUrl, tmpmask, res)

        for ide in filteredids('kegg', un, tmpmask) :
            qUrl = self.getLite + "search=" + \
                self.parent.urlSafe(ide) + "&searchCategory=DATABASE+LINK"
            fillresandmerge(qUrl, tmpmask, res)

        for ide in filteredids('inchi', un, tmpmask) :
            qUrl = self.getLite + "search=" + \
                self.parent.urlSafe(ide) + "&searchCategory=INCHI"
            fillresandmerge(qUrl, tmpmask, res)

        for ide in filteredids('smiles', un, tmpmask) :
            qUrl = self.getLite + "search=" + \
                self.parent.urlSafe(ide) + "&searchCategory=SMILES"
            fillresandmerge(qUrl, tmpmask, res)

        foundIds = res.keys()
        children = []
        for ch in foundIds :
            children = self.getChebiChildren(ch)
            for child in children :
                if not res.has_key(child):
                    res[child] = self.chebi2mask(mm, child)
        return(res)

    def chebi2mask (self, mm, chebiId) :
        """
        get a mask containing the info associated with a chebi id
        2. use getComplete to fetch the contents of the relevant entries
        """
        ba = chebiId
        un = mask({}, mm.idpatterns)
        qUrl = self.getComplete + "chebiId=" + str(chebiId)
        qRes = self.parent.getUrl(qUrl)
        if not qRes :
            return(mask({}))
        searchResults = xml.dom.minidom.parse(qRes)
        if not searchResults :
            return(mask({}))
        if searchResults.getElementsByTagName('ns1:return') :
            retList = searchResults.getElementsByTagName('ns1:return')[0]
        # found non-existent chebiId
        else :
            # this would have been reasonable but chebi has problems so that some
            # entries exist, but cant be fetched. valid chebi but we cant query
            # for it,
            # hence, skip completely
            # delmask = mask({})
            # delmask.append('chebi', chebiId)
            # mm.brandish(delmask)
            return(mask({}))
        # chebiid
        if retList.getElementsByTagName("ns1:chebiId") :
            chebiId = list(query.nodecontents(retList.getElementsByTagName("ns1:chebiId")))[0]
            un.append('chebi', chebiId.replace("CHEBI:", ""),
                      self.parent.confid, self.parent.sourceid)

        # smiles
        if retList.getElementsByTagName("ns1:smiles") :
            un.append('smiles', list(query.nodecontents(retList.getElementsByTagName("ns1:smiles")))[0],
                      self.parent.confid, self.parent.sourceid)

        # synonym
        if searchResults.getElementsByTagName('ns1:Synonyms') :
            syns = searchResults.getElementsByTagName('ns1:Synonyms')[0]
            for sy in list(query.nodecontents(syns.getElementsByTagName("ns1:data"))):
                un.append('synonym', sy, self.parent.mm.confidence['weak'], \
                              self.parent.sourceid)

        # inchi
        if retList.getElementsByTagName("ns1:inchi"):
            un.append('inchi', list(query.nodecontents(retList.getElementsByTagName("ns1:inchi")))[0],
                      self.parent.confid, self.parent.sourceid)

        # iupac
        if searchResults.getElementsByTagName('ns1:IupacNames'):
            syns = searchResults.getElementsByTagName('ns1:IupacNames')[0]
            for sy in list(query.nodecontents(syns.getElementsByTagName("ns1:data"))):
                un.append('iupac', sy, self.parent.confid,\
                              self.parent.sourceid)

        # kegg
        for ll in searchResults.getElementsByTagName('ns1:DatabaseLinks') :
            if list(query.nodecontents(ll.getElementsByTagName("ns1:type")))[0] == 'KEGG COMPOUND accession':
                for sy in list(query.nodecontents(ll.getElementsByTagName("ns1:data"))):
                    un.append('kegg', sy, self.parent.confid, self.parent.sourceid)

        # cas
        for ll in searchResults.getElementsByTagName('ns1:RegistryNumbers') :
            if list(query.nodecontents(ll.getElementsByTagName("ns1:type")))[0] == 'CAS Registry Number':
                for sy in list(query.nodecontents(ll.getElementsByTagName("ns1:data"))):
                    un.append('cas', sy, self.parent.confid, self.parent.sourceid)

        #pdb.set_trace()
        # formula
        if searchResults.getElementsByTagName('ns1:Formulae'):
            form = searchResults.getElementsByTagName('ns1:Formulae')[0]
            for sy in list(query.nodecontents(form.getElementsByTagName("ns1:data"))):
                un.append('formula', sy, self.parent.mm.confidence['weak'],\
                          self.parent.sourceid)

        return(un)

class parser :
    def __init__ (self, parent) :
        self.chebip = chebiParser(parent)
        parent.tables = self.chebip.queryTables
        self.parent = parent
        if not parent.boost:
            parent.boost = True
        parent.master = 'chebi'

    def process (self) :
        ll = True
        self.parent.mm.setTableWeak('formula')

        while ll :
            ll = self.parent.getLine()
            if not ll :
                break
            tmp = mask({})
            tmp.append('_id', ll)
            tmp2 = self.parent.mm.getMask(tmp)
            if not tmp2 :
                continue
            else :
                un = tmp2[0]
            try:
                chebiMasks = self.chebip.getChebiMasks(un, self.parent.mm)
            except:
                pdb.set_trace()
            if self.parent.mm.debug :
                print "#COMMENT mask " + str(ll) + " chebi " + str(chebiMasks.keys())
            chMask = mask({})
            for ch in chebiMasks.keys() :
                chebiMasks[ch].setAllAssoc(self.parent.mm.addAss())
                chMask.merge(chebiMasks[ch])
            self.parent.setMask(chMask, setass=False)
