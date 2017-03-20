import ConfigParser
import fileinput
import httplib
import os
import pdb
import re
import socket
import sys
import tempfile
import traceback
import urllib2
import xml.dom.minidom
from optparse import OptionParser

from SOAPpy import WSDL

import metmask
import metmask.dbi
import metmask.mask
import metmask.parse
from metmask.parse import main as pmain

reload(pmain)
mm = metmask.dbi.db(db="/tmp/balle", ask=False,\
                        debug=True,\
                        minoverlap=1)
myimporter = pmain.importer(mm, "chemspider", \
                                True,\
                                "cs",\
                                "good",\
                                token="c9c1cc45-e17e-4762-b5b9-cfd6d612ebc4")




csids = myimporter.parser.cs.SimpleSearch('glycine')
myimporter.parser.cs.GetCompoundInfo(mm, csids[0])

myimporter.parser.cs.ExtRefs(100, 'PubChem')



# from SOAPpy import SOAPProxy                      
# url = 'http://www.chemspider.com/Search.asmx'
# ns = "http://www.chemspider.com/"
# action = "http://www.chemspider.com/SimpleSearch"
# server = SOAPProxy(url, ns, action)
# server.config.dumpSOAPOut = 1 
# server.config.dumpSOAPIn = 1 
# server.SimpleSearch(query="asd", token="c9c1cc45-e17e-4762-b5b9-cfd6d612ebc4")



wsdl = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/eutils.wsdl'
server = WSDL.Proxy(wsdl)


bla = server.run_eSearch(tool='metmask', \
                             email='metmask-geek@lists.sourceforge.net', \
                             term="\"InChI=1/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h2,5,7-10H,1H2\"", usehistory='y', \
                             field="INCHI", retmax=10, \
                             db='pccompound') 


bla = server.run_eSearch(tool='metmask', \
                             email='metmask-geek@lists.sourceforge.net', \
                             term="c00072", usehistory='y', \
                             retmax=10, \
                             db='pccompound') 




bla = server.run_eSearch(tool='metmask', \
                             email='metmask-geek@lists.sourceforge.net', \
                             term="\"InChI=1S/H2O/h1H2\"", usehistory='y', \
                             field="INCHI", retmax=10, \
                             db='pccompound') 



myid = bla['IdList']['Id']

# returns a strange object? don't know how to navigate it..
res = server.run_eSummary(tool='metmask', \
                        email='metmask-geek@lists.sourceforge.net', \
                        retmax=10,\
                        db='pccompound',\
                        WebEnv=bla['WebEnv'],\
                        query_key=bla['QueryKey'])

# better
summaryUrl = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pccompound&tool=metmask&email=metmas-geek@lists.sourceforge.net&retmax=10&WebEnv=" + bla['WebEnv'] + "&query_key=" + bla['QueryKey']
summaryResults = xml.dom.minidom.parse(urllib2.urlopen(summaryUrl)).documentElement
