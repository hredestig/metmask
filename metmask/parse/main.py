from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import map
from builtins import str
from builtins import next
from builtins import object
import shlex
import socket
import six
from sqlite3 import ProgrammingError
if six.PY2:
    from urllib.request import urlopen as urlopen
    from  urllib.error import URLError as URLError
    import http.client as http_client
else:
    import urllib.request as urlopen
    import urllib.error as URLError
    import http.client as http_client

import metmask.parse
from metmask.parse import *

socket.setdefaulttimeout(10)


class parserError(Exception):
    """ raised when no suitable parser could be found or the parser had problems """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return (repr(self.value))


class fileFormatError(Exception):
    """ raised when the file to parse does not look like expected"""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return (repr(self.value))


def fixLine(ll, sep):
    """turn a delimited string in to a list
    with appropriate tokens
    eg '"a","2,2"' -> ['a', '2,2'] 

    Parameters:

    -`ll`, the string
    """
    # three exceptions to make any sure empty cells are not
    # missed, i.e. turn ',,' to ', ,'
    ll = ll.strip()
    ll = re.sub('^' + sep, ' ' + sep, ll)
    ll = re.sub(sep + '$', sep + ' ', ll)
    ll = re.sub(sep + sep, sep + ' ' + sep, ll)
    ll = re.sub(sep + sep, sep + ' ' + sep, ll)
    splitter = shlex.shlex(ll, posix=True)
    splitter.whitespace = sep
    splitter.whitespace_split = True
    res = list(splitter)
    # undo any damage we did
    res = [re.sub(sep + ' ' + \
                               sep, sep + \
                               sep, x) for x in res]
    return (res)


class importer(object):
    """A class for importing information to the metmask database
    """

    def __init__(self, mm,
                 parsertype, fileobj, source,
                 confidence='good', sep1=",", sep2="|", na="[nN][Aa]",
                 resolve=True, boost=False, master='unknown', token=None):
        if not re.match("_", parsertype):
            parsertype = "_" + parsertype

        if not parsertype in metmask.parse.PARSERS:
            raise parserError("Unknown parser")

        self.token = token
        """ the chemspider security token """
        self.master = master
        """ the master table for this imported source """
        self.tables = []
        self.sep1 = sep1
        """ primary separator (between tables) 123-23-3,c00001"""
        self.sep2 = sep2
        """ secondary separator (between indentifiers) 
        123-23-3|234-34-4"""
        self.na = na
        """ regexp to interprete as missing value """
        self.lineNum = 0
        """ current number of read lines """
        self.nentries = 0
        """ current number of submitted masks"""
        self.mm = mm
        """ the database """
        self.resolve = resolve
        """ should conflicting masks be resolved
        or merged directly  """
        self.boost = boost
        """ should masks only be added if they carry 
        overlap with something already in the database  """
        self.confidence = confidence
        """string describing the confidence"""
        self.confid = mm.addConf(confidence)
        """the confidenceid as defined by the db"""

        self.fileobj = fileobj
        # if we didnt get anything (just True), assume we loop over
        # all masks in the database
        if self.fileobj is True:
            self.fileobj = iter(mm.getAllMmids())
        # if we only got a string, assume it was a filename
        if 'next' not in dir(self.fileobj):
            self.fileobj = open(self.fileobj, 'r')
        """eventually an iterator giving the input"""

        self.source = source
        """string describing the source"""
        # ugly dynamic trick to get plug-in-esque parsers
        # the main action now happens with in 
        # importer.<_parser>.parser.process
        st = "self.parser = metmask.parse." + parsertype + ".parser(self)"
        exec(st, locals(), globals())

        self.sourceid = mm.addSource(source, master=self.master)
        """the sourceid as defined by the db"""

        # make sure that the necessary tables are in the db
        list(map(self.mm.createIdTable, self.tables))
        if master and master != 'unknown':
            self.mm.createIdTable(self.master)
        self.mm.connection.commit()

    def __del__(self):
        if 'close' in dir(self.fileobj):
            self.fileobj.close()
        try:
            self.mm.connection.commit()
        except ProgrammingError:
            pass

    def getLine(self, comment=None):
        """safe way to get a new line"""
        try:
            ll = next(self.fileobj)
            self.lineNum = self.lineNum + 1
            if comment:
                while str(ll).startswith(comment):
                    ll = next(self.fileobj)
                    self.lineNum = self.lineNum + 1
            if not re.match('\s', str(ll)):
                return (str(ll).strip())
            return (str(ll))
        except StopIteration:
            return ('')

    def setMask(self, ma, setass=True):
        """ set mask (ma) considering the settings of this parser
        Parameters:
        -`ma`: the mask to set
        -`setass`: should all associations be set to a new association code
        """
        goAhead = True
        if self.boost:
            mmids = self.mm.getMmid(ma)
            if not mmids:
                goAhead = False
        if goAhead:
            if setass:
                ma.setAllAssoc(self.mm.addAss())
            self.nentries = self.nentries + self.mm.setMask(ma, self.resolve)

    def urlSafe(self, string):
        """make string safe(r) for use as a url
        """
        string = string.replace("%", "%25")
        string = string.replace(" ", "+")
        string = string.replace("/", "%2F")
        string = string.replace("@", "%40")
        string = string.replace("=", "%3D")
        string = string.replace("[", "%5B")
        string = string.replace("]", "%5D")
        string = string.replace("+", "%2B")
        string = string.replace(":", "%3A")
        string = string.replace(";", "%3B")
        string = string.replace(",", "%2C")
        string = string.replace("&", "%26")
        string = string.replace("?", "%3F")
        string = string.replace("<", "%3C")
        string = string.replace(">", "%3E")
        string = string.replace("#", "%23")
        return (string)

    def getUrl(self, url):
        """get contents from url but do safely timeout if no response and
        ignore junk response
        """
        try:
            return urlopen(url)
        except http_client.BadStatusLine as inst:
            if self.mm.debug:
                print("#COMMENT bad response skipping")
                return (None)
        except URLError.URLError as inst:
            if self.mm.debug:
                print("#COMMENT no response skipping")
                return (None)
        except socket.timeout as inst:
            if self.mm.debug:
                print("#COMMENT no response skipping")
                return (None)
