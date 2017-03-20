""" 
A module to query online databases

Query tools
===========

a module to query other databases than the local metmask db. 

basic usage is if you have a mask instance but you lack some
identifiers in it, the you can use `boost` to add more identfiers to
it from a pubchem match (if one can be found of course). 

e.g::

  mm = dbi.db("yourdatabase")
  a = mm.simpleQuery('threitol, D-', what='synonym', to='mask')[0]
  query.fetch(mm, a)

"""

import pdb
import re
import urllib2
import xml.dom.minidom

from SOAPpy import WSDL

from mask import guessTable, mask

KNOWNTABLES = ['kegg', 'cas', 'cid', 'chebi', 'inchi', 'synonym']

MAX_SYN=10
MAX_RET=10

def nodecontents (nodes) :
    """generator for getting the contents of an xml node 

    Parameters :
    -`nodes`: a minidom node
    """
    for n in nodes :
        for val in n.childNodes :
            yield val.nodeValue


def pubchem2mask (mm, docSum, confidence='external') :
    """parse a minidom docSum node from pubchem and convert its contents
    to a mask, return the mask. cas, chebi and knapsack are good but
    pubchem does not indicate which indentifier is kegg so not
    everything that looks like a kegg identifier will be interepreted
    as one.

    Parameters :
    -`mm` : an instance of `dbi.db`
    -`docSum`: a minidom docSum node from pubchem
    -`confidence`: the desired confidence code
    """
    un = mask({}, mm.idpatterns)
    cid = nodecontents(docSum.getElementsByTagName("Id")).next()
    if mm.debug :
        print 'cid:' + cid,
    sourceid = mm.addSource('PubChem')
    confid = mm.addConf(confidence)
    un.append('cid', cid, confid, sourceid)

    for item in docSum.getElementsByTagName('Item'):
        if item.getAttribute('Name') == 'SynonymList' :
            synonyms = list(nodecontents(item.childNodes))
            for s in synonyms :
                # for clearing out links to kegg
                p = re.compile('(<a href[^>]*>)|(</a>)|(ligand)|(,)')
                s = p.sub('', s.lower()).strip()
                if 'kegg' in guessTable(s, mm.idpatterns) :
                    un.append('kegg', s, confid, sourceid)
                    continue
                if 'cas' in guessTable(s,  mm.idpatterns) :
                    un.append('cas', s, confid, sourceid)
                    continue
                # inchi identifiers from pccompound are corrupt
                # if re.match('inchi=', s) :
                #    un.append('inchi', s, confid, sourceid)
                # continue
                if re.match('chebi:', s) :
                    s = re.findall("chebi:(\d*)", s)[0]
                    un.append('chebi', s, confid, sourceid)
                    continue
                else :
                    un.append('synonym', s, mm.confidence['weak'], sourceid)
        if item.getAttribute('Name') == 'CommentList' :
            comments = list(nodecontents(item.childNodes))
            for c in comments :
                c = c.lower()
                x = re.findall("cas: (.*)", c)
                if x :
                    un.append('cas', x[0], confid, sourceid)
                x = re.findall("chebi: (.*)", c)
                if x :
                    un.append('chebi', x[0], confid, sourceid)
                x = re.findall("knapsack: (.*)", c)
                if x :
                    un.append('knapsack', x[0], confid, sourceid)
    return(un)

def kegg2mask (mm, entry, confidence='external') :
    """parse a kegg entry convert its contents to a mask, return the
    mask. 

    Parameters :
    -`mm` : an instance of `dbi.db`
    -`entry`: string as gotten from KEGGs bget
    -`confidence`: the desired confidence code
    """
    un = mask({},  mm.idpatterns)

    keggId = re.findall('ENTRY +([C|D]\d{5})', entry)
    if not keggId :
        raise NameError, 'weird looking kegg entry'
    keggId = keggId[0]

    sourceid = mm.addSource('KEGG')
    confid = mm.addConf(confidence)
    un.append('kegg', keggId, confid, sourceid)
    
    entsplit = entry.split('\n')
    entsplit.reverse()
    e = entsplit.pop()
    while e :
        # search in name field for synonyms
        if re.match('NAME', e) :
            syn = re.findall('NAME +([^;]+)', e)[0]
            un.append('synonym', syn, confid, sourceid)
            e = entsplit.pop()
            while re.match(' +.+', e) :
                syn = re.findall(' +([^;]+)', e)[0]
                un.append('synonym', syn, confid, sourceid)
                e = entsplit.pop()

        # search in dblink field for other identifiers
        if re.match('DBLINKS', e) :
            # cas
            x = re.findall('CAS: (\d{2,7}-\d{2}-\d{1})', e)
            if x :
                un.append('cas', x[0], confid, sourceid)
            # cid
            x = re.findall('PubChem: (\d+)', e)
            if x :
                un.append('cid', x[0], confid, sourceid)
            # chebi
            x = re.findall('ChEBI: (\d+)', e)
            if x :
                un.append('chebi', x[0], confid, sourceid)
            while re.match(' +.+', e) :
                # cas
                x = re.findall('CAS: (\d{2,7}-\d{2}-\d{1})', e)
                if x :
                    un.append('cas', x[0], confid, sourceid)
                # cid
                x = re.findall('PubChem: (\d+)', e)
                if x :
                    un.append('cid', x[0], confid, sourceid)
                # chebi
                x = re.findall('ChEBI: (\d+)', e)
                if x :
                    un.append('chebi', x[0], confid, sourceid)
                e = entsplit.pop()
        e = entsplit.pop()
            
    
    return(un)

def getPubchem (mm, un) :
    """get a bunch of pubchem based masks for the identifiers given in the
    current mask. Return the masks in a list if it managed to add
    something to the mask, [] otherwise

    Parameters:
    
    -`mm`: an instance of `dbi.db`
    -`un`: an instance of `mask`
    """
    # tables other than synonym that are known for pubchem
    tables = ['kegg', 'cas', 'cid', 'chebi', 'inchi']

    res = []
    if not any(map(lambda x: x in KNOWNTABLES, un.getTables())) :
        return(res)
    esUrl = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&retmax=" + str(MAX_RET) + "&tool=metmask&email=metmask-geek@lists.sourceforge.net&usehistory=y&term="

    # add synonym to specific field and the rest just as 'term'
    if un.hasTable('synonym') :
        identifiers = un.getIdentifiers('synonym')
        for i in range(0, min(MAX_SYN, len(identifiers))) :
            identifiers[i] = re.sub('\s', '+', identifiers[i])
            esUrl = esUrl + identifiers[i] + "[synonym]+OR+" 

    for tab in tables :
        if un.hasTable(tab) :
            for ide in un.getIdentifiers(tab) :
                ide = re.sub('\s', '+', ide)
                esUrl = esUrl + ide + "[All Fields]+OR+"
    esUrl = re.sub("\+OR\+$", "", esUrl)
    searchResults = xml.dom.minidom.parse(urllib2.urlopen(esUrl))
    
    # fetch the results
    webenv = list(nodecontents(searchResults.getElementsByTagName('WebEnv')))[0]
    querykey = list(nodecontents(searchResults.getElementsByTagName('QueryKey')))[0]
    summaryUrl = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pccompound&tool=metmask&email=metmas-geek@lists.sourceforge.net&retmax=" + str(MAX_RET) + "&WebEnv=" + webenv + "&query_key=" + querykey
    summaryResults = xml.dom.minidom.parse(urllib2.urlopen(summaryUrl)).documentElement
    
    # guess on the first compatible mask, if there is one, return it
    candidates = summaryResults.getElementsByTagName('DocSum')
    for c in candidates :
        pcmask = pubchem2mask(mm, c)
        res.append(pcmask)
    return(res)

def getKegg(mm, un) :
    """get a bunch of KEGG based masks for the identifiers given in the
    current mask. Return the masks in a list if it managed to add
    something to the mask, [] otherwise

    Parameters:
    
    -`mm`: an instance of `dbi.db`
    -`un`: an instance of `mask`
    """
    res = []
    if not any(map(lambda x: x in KNOWNTABLES, un.getTables())) :
        return(res)

    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    kegg = WSDL.Proxy(wsdl)
    
    
    keggIds = []
    if un.hasTable('synonym') :
        for s in un.getIdentifiers('synonym') :
            keggIds.extend(kegg.search_compounds_by_name())
    if len(keggIds) > MAX_RET:
        keggIds = keggIds[0:(MAX_RET - 1)]
        
    if un.hasTable('kegg'):
        keggIds = keggIds + un.getIdentifiers('kegg')

    entries = []
    for k in keggIds :
        if not re.match('cpd:', k):
            hit = kegg.bget('cpd:' + k)
        else :
            hit = kegg.bget(k)

        if hit :
            entries.append(hit)

    res = []
    for e in entries :
        keggmask = kegg2mask(mm, e)
        res.append(keggmask)

    return(res)
        
def fetch (mm, un, internal, to) :
    """ try to merge the input mask with a new information from the
    web. Return the new mask (which is compatible with existing
    information) if successful, False otherwise. This is the main
    function to use from outside the query module.

    Parameters:
    
    -`mm`: an instance of `dbi.db`
    -`un`: an instance of `mask`
    -`internal`: logical, is the query mask a mask that was already in
     the database, in which we need the result to be compatible with
     that mask?  If it is not internal we do not check compatibility
     (not necessary as there won't be any conflicts anyway)
    """
    if to[0] == 'mask':
        to = KNOWNTABLES

    masks = []
    if any(map(lambda x:x in 'kegg', to)) :
        masks.extend(getKegg(mm, un))

    if any(map(lambda x:x in ['synonym', 'chebi', 'cid', 'cas', 'inchi'], to)) :
        masks.extend(getPubchem(mm, un))

    # TODO
    # this should somehow be done better so that obtained masks are
    # 'clustered' themselves before being checked with the database

    for pcmask in masks :
        suggestion = un.resolve(pcmask, mm)
        # check if we are allowed to merge and if, then do so
        # other wise don't
        if suggestion == 'merge' or not internal:
            return(pcmask)
    # means that we didn't merge
    return(False)
