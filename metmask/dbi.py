""" 
A module to manage all transcations with the metmask database.

How to use this module
======================

1. Get a connection to your database (just a file location or possible
   ':memory:' if you do intend to keep the information
   
   ``mm = dbi.db('/home/user/metmask.db')``

2. Insert with `setMask` (see the `mask` module) 

3. Query information with the `simpleQuery` interface

4. Close the database either explicitly with `close` or let python
   take care of it with just `del mm`

Database tables
===============

- version, the current version number to make sure we do not try
  to query an incompatible database
- confidence, holds the confidence codes 
- sources, holds the sources codes 
- assidx, just an column index (could just as well have been
  part of the associations table
- assocaition, holds information about all associations, their
  confidence and where they came from
- idtables, a table specifiying which identifier tables the
  database currently posses
- _id, basically also just an identifier table but this holds
  the master mask identifier
- the rest are tables that holds the identifiers and are named
  exactly like the identifier it self

Copyright (C) Henning Redestig
2009
See COPYING.txt for licensing details.

"""

import pdb, re, sys, os
import sqlite3
from mask import mask
import query
import metmask

def determine_path ():
   """Borrowed from wxglade.py"""
   try:
      root = __file__
      if os.path.islink (root):
         root = os.path.realpath (root)
      return os.path.dirname (os.path.abspath (root))
   except:
      print "I'm sorry, but something is wrong."
      print "There is no __file__ variable. Please contact the author."
      sys.exit ()

class recursionError(Exception) :
   """this error indicates that when setting a new mask, the program got
   stuck in a recursion, attempting to merge too many mask. This
   should not happen unless the database is corrupt or the user tries
   to set a mask that holds too many common identifiers.
   """
   def __init__(self, value):
       self.value = value
   def __str__(self):
       return(repr(self.value))

class dbError(Exception) :
   """generic database error that wasn't due to syntax problems
   """
   def __init__(self, value):
       self.value = value
   def __str__(self):
       return(repr(self.value))

class db :
   """ The object that should do all transcation to the database. The use
   of a single db object for each database ensures that no db-locks
   should occur and that connections are not opened and closed too often.
   """
   def __init__(self, \
                   db, \
                   ask=False,\
                   debug=True,\
                   minoverlap=2) :
      """ Initialize a `db` object.
      
      Parameters
      ----------
      
      -`db`, string indicating where the database is or should be
       stored. Use ':memory:' to keep everything in the working memory.

      - `ask`, logical indicating if the database should expect
        feedback from the user or not when deciding how to merge
        masks.

      - `debug`, logical for outputting more information.

      """
      
      self.db = db
      self.ask = ask
      self.debug = debug

      self.ifconflict = 'm'
      """ `string` indicating the default action when merging (merge, prune) """
      self.confidence = {}
      """ `dictionary` holding the current known confidence codes and the
      indices in the database """
      self.idpatterns = {'kegg'      : '[cCdD]\d{5}$',
                         'cas'       : '\d{1,7}-\d{2}-\d{1}$',
                         '_id'       : '\d+$',
                         'synonym'   : '.+',
                         'mpimp'     : '.+',
                         'riken'     : '.+',
                         'cycdb'     : '.+',                
                         'rlib'      : '[a-z|A-Z]{0,1}\d{4}',
                         'kappav'    : '[kK][pP][cC]\d{5}',
                         'knapsack'  : '[cC]\d{8}',
                         'tsuruoka'  : '[tT]-\d+$',
                         'cqid'      : 'cq_\d{5}$',
                         'cid'       : '\d+$',
                         'chebi'     : '(CHEBI:){0,1}\d+$',
                         'default'   : 'synonym',
                         'inchi'     : 'InChI=.+',
                         'metlin'    : '\d+$',
                         'sid'       : '\d+$',
                         'hmdb'      : 'hmdb\d+$',
                         }
      """ dictionary which holds regular expression that 
      all identifiers for the respective table must match """
      self.METMASK_CONF = 0
      """ reserved confidence index for internal associations """
      self.FORBIDDENTABLENAMES = ['ALL', 'bioconductor', 'all','mask','graph', 'unknown']
      """ reserved strings """
      self.METMASK = 'metmask'
      self.MAIL = 'metmask-geek@lists.sourceforge.net'
      self.METMASKID  = 0
      """ reserved source string  """
      self.MIN_OVERLAP = int(minoverlap)
      """`int` the mininum number of tables this mask must match with other
      masks in order to merged directly """
      self.depth = 0
      """ the current level of recursion, used to break a transcation before
      recurses too deep """
      self.rw = False
      self.c = None
      """ indicates if we have rw access or not """
      if not os.access(self.db, os.F_OK) :
         try: 
            os.mknod(self.db)
         except :
            try:
               open(self.db, 'a').close()
            except :
               raise dbError, 'can not create db @ ' + self.db
      if not os.access(self.db, os.R_OK):
         raise dbError, 'can not access db @ ' + self.db
      self.rw = os.access(self.db, os.W_OK) 
      self.connection = sqlite3.connect(self.db)
      """ `sqlite3.connection` the current connection to the database """
      self.connection.text_factory = str
      """ see sqlite3 documentation """
      self.c = self.connection.cursor()
      """ `sqlite3.cursor` the object that can write the db """
      if self.debug:
         print "#COMMENT Welcome to metmask v" + metmask.__version__
      if not self.getIdTables() :
         self.setup() 
         if self.debug:
            print "started new database"
      else:
         # check that we have a legal version (same middle version
         # number), other version numbers indicate database
         # structural changes
         self.c.execute("SELECT number FROM version")
         version = self.c.fetchone()[0]
         if metmask.__version__.split('.')[1] != version.split('.')[1]:
            raise dbError, 'existing database was created with version ' + \
                str(version) + \
                ', the database format has changed, re-create the database or use a matching metmask version'

      self.updateConfidence()
      self.updateIdpatterns()
      
   def __del__(self) :
      """ commit and close and connection """
      try :
          self.close()
      except sqlite3.ProgrammingError, inst :
          if inst[0] != 'Cannot operate on a closed database.' :
              raise inst

   def updateIdpatterns (self) :
      """ update the patterns dictionary """
      self.c.execute("SELECT * FROM idpatterns")
      tmp = self.c.fetchall()
      for t in tmp :
         self.idpatterns[str(t[0])] = t[1]
      if self.rw :
         self.c.execute("DELETE FROM idpatterns")
         for k in self.idpatterns.keys() :
            self.c.execute("INSERT INTO idpatterns VALUES (?,?)", (k, self.idpatterns[k]))

   def updateConfidence (self) :
      """ update the confidence dictionary with indices and the codes """
      self.c.execute("SELECT * FROM confidence")
      for co in self.c.fetchall() :
         self.confidence[str(co[1])] = int(co[0])

   def close(self):
      """ commit and close the connection """
      if self.c :
         self.updateForBioC()
         self.connection.commit()
         self.connection.close()
         if self.debug:
            print "#COMMENT closed metmask database"
      
   def setup (self):
      """ Create the necessary tables for and insert basic information that
      must be in all metmask databases.
      """
      try:
         # misc tables
         self.c.execute("CREATE TABLE version (number TEXT)")
         self.c.execute("INSERT INTO version VALUES (:num)", {'num':metmask.__version__})
         self.c.execute("CREATE TABLE metadata (name TEXT, value TEXT)")
         self.c.execute("CREATE TABLE map_counts (map_name TEXT, count INTEGER)")
         self.c.execute("CREATE TABLE  idpatterns (name TEXT, pattern TEXT)")
         self.c.execute("CREATE TABLE assidx (assid INTEGER PRIMARY KEY)")

         # confidence
         st  = "CREATE TABLE confidence ( confid INTEGER PRIMARY KEY, description TEXT UNIQUE )"
         self.c.execute(st)
         self.c.execute("INSERT INTO confidence VALUES (?, ?)", (self.METMASK_CONF, 'metmask-internal'))
         self.c.execute("INSERT INTO confidence VALUES (?, ?)", (metmask.NEVERMERGE_CONF, 'nevermerge'))
         self.c.execute("INSERT INTO confidence VALUES (?, ?)", (metmask.WEAK_CONF, 'weak'))

         # sources
         self.c.execute("CREATE TABLE sources ( sourceid INTEGER PRIMARY KEY, name TEXT, master TEXT REFERENCES idtables (name))")
         self.c.execute("INSERT INTO sources VALUES (?, ?, ?)", (self.METMASKID, self.METMASK, "_id"))

         # standards id tables
         self.c.execute("CREATE TABLE _id (_id INTEGER PRIMARY KEY, mmid INTEGER, conf INTEGER, src INTEGER, qc INTEGER DEFAULT 0, ass INTEGER DEFAULT 0, UNIQUE (_id, qc))") # mmid is duplicate of _id
         self.c.execute("CREATE TABLE preferred (_id INTEGER UNIQUE REFERENCES _id (_id), preferred TEXT UNIQUE, conf INTEGER, src INTEGER, qc INTEGER, ass INTEGER DEFAULT 0)")

         # idtables
         self.c.execute("CREATE TABLE idtables ( tableid INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, weak INTEGER DEFAULT 0 )")
         self.c.execute("INSERT INTO idtables VALUES (0, \"_id\", 0)")
         self.c.execute("INSERT INTO idtables VALUES (1, \"preferred\", 0)")
      except sqlite3.OperationalError :
          raise dbError, "Something went wrong, is the db already initialized?"

   def updateForBioC (self) :
      """ update map_counts and etc so bioconductor is happy
      """
      if not self.rw :
         return 0
      self.c.execute("DELETE FROM map_counts")
      tabs = self.getIdTables()
      total = 0
      for t in tabs:
         self.c.execute("SELECT count(*) FROM " + t)
         num = self.c.fetchall()[0][0]
         total = total + num
         self.c.execute("INSERT INTO map_counts VALUES (?,?)", (str(t).upper(), int(num)))
      self.c.execute("INSERT INTO map_counts VALUES (?,?)", ("TOTAL", int(total)))

      self.c.execute("DELETE FROM metadata")
      self.c.execute("INSERT INTO metadata VALUES (?, ?)", (str("DBSCHEMAVERSION"), str("2.0")))
      self.c.execute("INSERT INTO metadata VALUES (?, ?)", (str("DBSCHEMA"), str("METMASK_DB")))
      self.c.execute("INSERT INTO metadata VALUES (?, ?)", ("METMASK_VERSION", str(metmask.__version__)))
      for k in self.confidence.keys():
         self.c.execute("INSERT INTO metadata VALUES (?, ?)", ("CONFIDENCE_" + str(k).upper(), self.confidence[k]))
      self.c.execute("SELECT * FROM sources")
      tmp = self.c.fetchall()
      for t in tmp :
         self.c.execute("INSERT INTO metadata VALUES (?, ?)", ("SOURCE_"+t[1].upper(), t[0]))
      strong = self.getIdTables(weak=False)
      for t in tabs:
         weak = "WEAK"
         if t in strong:
            weak = "STRONG"
         self.c.execute("INSERT INTO metadata VALUES (?, ?)", (t.upper(), weak))

   def createIdTable (self, name, weak=0) :
      """ Create a new identifier table. No return value.
      
      Parameters:
      - `name`: string specifying the name of the identifier.
      """
      name = name.rstrip()
      if re.match(".*\..*", name) :
         raise dbError, "dot in table name " + name + ". Specify other name"
      if name in self.FORBIDDENTABLENAMES :
          raise dbError, "name of table must not be one of" + str(self.FORBIDDENTABLENAMES)
      stTable = "CREATE TABLE " + name + " ( _id INTEGER REFERENCES _id (_id), " + \
          name + \
          " TEXT, conf INTEGER, src INTEGER, qc INTEGER, ass INTEGER DEFAULT 0, UNIQUE (_id, ass, " + name + "), UNIQUE (" + name + ", qc))"
      try:
         self.c.execute(stTable)
         self.c.execute("CREATE INDEX " + name + "_idx" + " ON " + name + " ( " + name + " )")
         self.c.execute("INSERT INTO idtables VALUES (?, ?, ?)", (None, name, weak))
      except sqlite3.OperationalError, inst :
         if not re.match(".*already exists$", str(inst)) :
             raise inst

   def setTableWeak(self, name) :
      """ grade an id-table as weak (only relevant for bioconductor)
      
      Parameters:
      -`name` : the name of the table
      """
      if not name in self.getIdTables() :
         raise dbError, "trying to weaken a non-existent table: " + str(name)
      self.c.execute("UPDATE idtables SET weak = 1 WHERE name = :name", {'name':name})

   def addAss(self) :
      """ add an association to the database and get the new association
      index back
      """
      st = "INSERT INTO assidx VALUES (:ass)"
      self.c.execute(st, {'ass':None})
      return(self.c.lastrowid)

   def addMm(self) :
      """ add a new mask and get it the index back 
      """
      self.c.execute("INSERT INTO _id VALUES (?, ?, ?, ?, ?, ?)",(None, 0, self.METMASK_CONF, self.METMASKID, 0, 0))
      _id = self.c.lastrowid
      self.c.execute("UPDATE _id SET mmid = :_id WHERE _id = :_id", {'_id':_id})
      return(_id)

   def dropPreferred(self, _id) :
      """ delete any preferred labels that might be associated with this
      mask

      Parameters:
      - `_id`: the _id 
      """
      ma = mask({})
      ma.append('_id', _id)
      existing = self.getMask(ma)
      if existing:
         if existing[0].hasTable('preferred') :
            prefs = existing[0].getIdentifers('preferred')
            for p in prefs :
               self.c.execute('DELETE FROM preferred WHERE _id = :_id', {'_id':_id})
            
   def addConf (self, name) :
      """  add or fetch a entry from/to the confidence table, get the
      created/existing confidence index back 
      """
      if self.confidence.has_key(name) :
         return(self.confidence[name])
      try :
         self.c.execute("INSERT INTO confidence VALUES (?, ?)", (None, name))
      except sqlite3.IntegrityError, inst :
         raise inst
      self.updateConfidence()
      return(self.addConf(name))

   def addSource (self, name, add=True, master='unknown') :
      """  add or fetch a entry from/to the sources table, get the
      created/existing source index back 
      """
      if name in self.getIdTables() :
          raise dbError, "source name must not equal a table name"
      self.c.execute("SELECT sourceid FROM sources WHERE name = :name", {'name':name})
      
      sourceid = self.c.fetchall()
      if sourceid or not add :
          return(sourceid[0][0])
      else :
          self.c.execute("INSERT INTO sources VALUES (?, ?, ?)", (None, name, master))
          return(self.c.lastrowid)

   def sourceid2master(self, seq) :
      """ convert a list of source id's to source names 

      Parameters:
      - `seq`, a list of source names
      """
      res = []
      st = "SELECT master FROM sources WHERE sourceid = ? "
      for s in seq :
         self.c.execute(st, (s,))
         tmp = self.c.fetchall()
         for t in tmp :
            res = res + [t[0]]
      return(res)

   def sourceid2source(self, seq) :
      """ convert a list of source id's to source names 

      Parameters:
      - `seq`, a list of source names
      """
      res = []
      st = "SELECT name FROM sources WHERE sourceid = ? "
      for s in seq :
         self.c.execute(st, (s,))
         tmp = self.c.fetchall()
         for t in tmp :
            res = res + [t[0]]
      return(res)

   def getWeakTables(self) :
      alltab = self.getIdTables(True)
      strong = self.getIdTables(False)
      return(filter(lambda x : x not in strong, alltab))

   def getIdTables(self, weak=True) :
      """  get a list with the names of the tables that hold identifiers
      
      Parameters:
      -`weak`: fetch weak tables as well
      """
      try :
         if not weak:
            self.c.execute("SELECT name FROM idtables WHERE weak != 1")
         else:
            self.c.execute("SELECT name FROM idtables")
         tables = map(lambda x:x[0], self.c.fetchall())
         return(tables)
      except sqlite3.OperationalError, inst:
         if re.match('no such table', str(inst)) :
            return(None)
         else :
            raise inst

   def getMmid (self, un, weak=False, wildcards=False) :
      """ get the corresponding mmids pointing to anything in this
      mask. Return a vector of mmids.
      
      Parameters:
      - `un` : a mask
      """

      tables = un.getTables()
      res = set()
      for tab in tables :
         st = "SELECT _id FROM " + tab + " WHERE " + tab + "= :identifier"
         if wildcards  :
            st = "SELECT _id FROM " + tab + " WHERE " + tab + " LIKE :identifier"
         if not weak:
            st = st + " AND conf !=" + str(metmask.WEAK_CONF)
         identifiers = un.getIdentifiers(tab, weak=weak)
         for ide in identifiers :
            try :
               self.c.execute(st, {'identifier':ide}) 
               tmp = self.c.fetchall()
               if tmp : 
                  for t in tmp :
                     res.add(t[0])
            except sqlite3.OperationalError, inst :
               if re.match("no such table", str(inst)) :
                  pass
               else :
                  raise
      # loose the set feature
      res = map(lambda x: x, res)
      return(res)

   def getMask (self, un, weak=False, wildcards=False) :
      """ get all masks that point to anything in the given mask. 

      Parameters:
      - `un`, an instance of `mask.mask`
      """
      mmids = self.getMmid(un, weak=weak, wildcards=wildcards)
      result = []
      tabs = self.getIdTables()
      for mm in mmids :
         newmask = mask({})
         for ta in tabs:
             st = "SELECT * FROM " + ta + " WHERE _id = " + str(mm)
             if not weak:
                st = st + " AND conf != " + str(metmask.WEAK_CONF)
             try :
                self.c.execute(st)
                tmp = self.c.fetchall()
             except sqlite3.OperationalError, inst :
                if re.match('no such table', inst) :
                   pass
                else:
                   raise inst
             if tmp :
                for t in tmp :
                   ide = t[1]
                   if ta == '_id' :
                      ide = int(ide)
                   newmask.append(ta, ide, t[2], t[3], t[5])
         result.append(newmask)
      return(result)

   def insertToMask(self, table, _id, value, conf, src, retry=True, nozero=False, ass=0) :
      """ Insert something to an id-table
      
      Parameters:
      -`table`: the table to insert to
      -`_id` : the _id (aka mmid)
      -`value` : the value to associate _id with (eg cas number)
      -`conf` : the confidence code
      -`src` : the source code
      """

      if table not in self.getIdTables(weak=True):
         raise dbError, "trying to insert to unknown table: " + str(table)
      if table in self.getWeakTables() :
         conf = metmask.WEAK_CONF
      if table == '_id' or table == 'preferred':
         conf = self.METMASK_CONF
      if conf == metmask.WEAK_CONF or nozero:
         self.c.execute("SELECT count(*) FROM " + table)
         qc = self.c.fetchall()[0][0] + 1
      else:
         qc = 0
      st = "INSERT INTO " + table + " VALUES (:id, :value, :conf, :src, :qc, :ass)"
      try :
         self.c.execute(st, {'id':_id, 'value':value, 'conf':conf, 'src':src, 'qc':qc, 'ass':ass})
      except sqlite3.IntegrityError, inst :
         if re.match("columns.+are not unique", str(inst)):
            if retry:
               # the ugly hack to allow strong identifier to come from
               # several sources
               stx = "SELECT DISTINCT _id FROM " + table + " WHERE " + table + " = \"" + value + "\" AND conf !=" + str(metmask.WEAK_CONF)
               self.c.execute(stx)
               checkOk = map(lambda x:x[0] == _id, self.c.fetchall())
               if all(checkOk) :
                  self.insertToMask(table, _id, value, conf, src, False, True, ass=ass)
            else:
               pass

   def setMask (self, un, resolve=True) :
      """ Try to input all information in this mask to the database. The db
      is first queried for any conflicting masks and the new
      information is then either merged to an existing mask or
      conflicting information is pruned. See manual.txt for further
      clarification of how this is decided.
      
      Parameters:

      - `un`: a instance of `mask.mask`
      - `resolve`: should conflicts try to be resolved or just over-run with 'merge'
      """
      # some globals for this function
      mmidsToMerge = []
      def panicCheck() :
         self.depth = self.depth + 1
         if self.depth > 100 :
            return(True)
         else:
            return(False)
         
      # if you add something, conditions are changed and we might get
      # stuck in a loop where synonym is removed due to conflicet then
      # re-added here and...  unless we only do it once... HACK!
      # preferred must also be a strong synonym
      if self.depth == 0:
         if un.hasTable('preferred') :
            mpref = un.getIdentifiers('preferred')[0]
            there = False
            if un.hasTable('synonym'):
               if mpref in un.getIdentifiers('synonym', weak=False):
                  there = True
            if not there:
               un.append('synonym', mpref, un.getConfidence('preferred', mpref)[0], \
                            un.getSource('preferred', mpref)[0])

      # ignore pointless masks, (only maps to itself)
      #       un.show()
      #       print "-------------------"
      #if un.nid() < 2 :
      if un.isEmpty() :
         self.depth = 0                           # no risk to get stuck agaion
         if self.debug:
            print "#COMMENT Ignoring pointless mask"
         return(0)

      # see if we already know something about its content, otherwise
      # we get []
      _id = self.getMmid(un)
      conflict = False
      if _id :
         conflict = True
   

      if conflict :
         # There are mmids that already point to identifiers behind
         # this mask, such identifiers are not allowed.  solution is
         # to either:
         #
         # - merge the new mask with the existing mask.
         #
         # - create a new mask but prune away the offending identifiers
         #   (what we do if there is an error in the input)
         #

         # get the first conflicting masks
         tmpmask = mask({})
         tmpmask.append('_id', _id[0])
         cnf = self.getMask(tmpmask, weak=True).pop() # get the full mask

         if panicCheck() :
            if self.debug :
               pdb.set_trace()
            else :
               raise recursionError, 'bogging down..'

         # try to figure out what to do with the rest:
         if resolve: 
            possibility = un.resolve(cnf, mm=self)
         else :
            possibility = 'merge'
         if not possibility :
            if self.debug :
               print "#COMMENT keeping existing mask..",
            un.subtract(cnf)                      # subtract non-weak ids
            self.setMask(un, resolve)
            return(0)
         else :
            self.ifconflict = possibility[0]
            if self.debug:
                print "#COMMENT " + str(possibility)
         
         # >>>>Get user decision if asked for
         if self.ask:
            print "***** conflicting mask *******"
            cnf.show(mm=self)
            print "***** mask to insert *********"
            un.show(mm=self)
            print "***** overlap ****************"
            intersection = un.intersect(cnf)
            intersection.show(mm=self)
            print "default:" + possibility
            choice =  raw_input("[m]erge/[p]rune and merge : ")
            if len(choice) > 0 :
               self.ifconflict = choice
         if self.ifconflict[-1] == '!' :
            self.ask = False
            self.ifconflict = self.ifconflict[0]
         # done asking>>>>>>>>>>>>>>>>>>>>>>>

         # merge the new information into the existing one
         if self.ifconflict in ['','m'] :
            if self.debug:
               print "#COMMENT merging to mask..",
            un.merge(cnf)                 # got an _id.. , merge weak as well
            self.dropMask(cnf)                    # drop weak as well
         if self.ifconflict == 'p' :
            # now prune away what is left
            un.subtract(cnf) # lost the _id.. .. but no overlap
         self.setMask(un, resolve)# ..anyway, re-submit
         return(1)

      elif not conflict :

         if un.hasTable('_id') :
            # make sure we remove all mmids but the first one, the
            # other were already removed from the db
            allMmids = un.getIdentifiers('_id')
            mmidsToMerge = un.getIdentifiers('_id')
            while len(allMmids) > 1:
               un.delIdentifiers('_id', allMmids.pop())
            
         # get a new mmmid
         if not un.hasTable('_id') :
            # no _id matches the given identifiers, we have found a new compound
            # add to the 'preferred' table
            _id = self.addMm()
            # adding a clause here to set conf to nevermerge if any
            # nevermerge, we could probably use the weakenIdentifier again
            if un.anyNeverMerge():
               un.append('_id', _id, 0, metmask.NEVERMERGE_CONF)
            else:
               un.append('_id', _id, 0, self.METMASK_CONF)
            self.fixPreferred(un)

      self.depth = 0
      # get the (old or new) _id without its dimension
      _id = un.getIdentifiers('_id')[0]
      # replace
      if self.debug :
         print "inserting mask:" + str(_id)

      # first drop any preferred we already have in the database, we
      # set it back again
      if un.hasTable('preferred') :
         self.fixPreferred(un)
         self.dropPreferred(_id)

      for tab in un.getTables() :
         for ide in un.getIdentifiers(tab) :
            try :
               conf = un.getConfidence(tab, ide)
               src = un.getSource(tab, ide)
               ass = un.getAssoc(tab, ide)

               if len(conf) != len(src) or len(conf) <= 0  or len(src) <= 0 :
                  raise dbError, 'confidence and source not defined or not of equal length'
               
               for i in range(0, len(conf)) :
                  self.insertToMask(tab, _id, str(ide), conf[i], src[i], ass=ass[i])

            except sqlite3.IntegrityError, inst :
               raise
      return(1)

   def fixPreferred(self, un) :
      """ make sure that the preferred string is present and not
      available in the preferred table.

      Parameters:
      -`un` : a mask
      """
      if not un.hasTable('_id') :
         raise dbError, "get an _id before fixing the preferred"
      if un.hasTable('preferred') :
         self.c.execute("SELECT _id FROM preferred WHERE preferred = :pref", \
                           {'pref':un.getIdentifiers('preferred')[0]})
         checkOk = map(lambda x:x[0] == un.getIdentifiers('_id')[0], self.c.fetchall())
         if all(checkOk) :
            return(1)
         else:
            un.delTable('preferred')
            self.fixPreferred(un)
      if not un.hasTable('preferred') :
         if un.hasTable('synonym') :
            if un.getIdentifiers('synonym', weak=False):
               pref = un.getIdentifiers('synonym', weak=False)[0]
            else:
               pref = un.getIdentifiers('synonym')[0]
            so = un.getSource('synonym', pref)[0]
            cn = un.getConfidence('synonym', pref)[0]
         else :
            pref = un.getIdentifiers(un.getTables()[0])[0]
            so = un.getSource(un.getTables()[0], pref)[0]
            cn = un.getConfidence(un.getTables()[0], pref)[0]
            #this was meant to it recognizable if pref was
            #auto-generated but it also conteracts the purpose of
            #pref as it the string becomes longer
            #pref = 'mm:' + pref
         pref = str(un.getIdentifiers('_id')[0]) + ':' + str(pref)             # make pref unique
         un.append('preferred', pref, cn, so)
         self.fixPreferred(un)


   def dropIdentifiers(self, un) :
      """ drop all identifiers in the current mask. (Not the same as
      dropMask which deletes a single whole mask)
      
      Parameters:
      -`un` : a mask object containign the identifiers to drop.
      """
      tables = un.getTables()
      for tab in tables :
         identifiers = un.getIdentifiers(tab)
         for ide in identifiers :
            st1 = "DELETE FROM " + tab + " WHERE " + tab + "  = :ide"
            try: 
               self.c.execute(st1, {'ide':str(ide)})
            except sqlite3.OperationalError, inst :
               raise inst

   def dropMask (self, un) :
      """ drop all information related to this mask. Warning: annotations
      are kept. Make sure they are updated after deletion if necessary.

      Parameters:
      - `un`: instance of mask.mask
      """
      mmids = self.getMmid(un)
      if not mmids :
         raise dbError, 'no such mask'
      tables = self.getIdTables()
      for mm in mmids :
         for tab in tables :
            st1 = "DELETE FROM " + tab + " WHERE _id = :mm"
            try: 
               self.c.execute(st1, {'mm':str(mm)})
            except sqlite3.OperationalError, inst :
               raise inst

   def simpleQuery(self, name, what='_id', to='_id',\
                      external=True, \
                      learn=True,\
                      weak=True, outmask=False, wildcards=False): 
      """ Query the database for a single or set of identifiers.  Return a
      vector of masks. Or the empty vector if not hits were found.

      Parameters : 
      
      - `name`, the identifier(s) eiher as a string or list of strings.
      - `what`, the type of the identifier, if a list each identifier
      is interpreted as potentially one of each 
      - `to`, string indicating the desired return type. Either return a mask or the
      desired table
      - `external`, logical, should try internet if the database
      didn't have the desired identifier
      - `weak`, logical, should we return weak identifiers as well
      - `outmask`, logical, should we return a mask or formatted list
      """
      if not weak:
         weak = False
      originalTo = to
      un = mask({})

      # make sure that inputs are proper lists
      if to[0] == 'ALL':
         to = self.getIdTables(weak)
      if re.match("~", to[0]):
         to = filter(lambda x:x not in map(lambda x:x.replace("~", ""), to),\
                        self.getIdTables(weak))
      if not isinstance(to, list) :
         to = [to]
      if not isinstance(what, list) :
         what = [what]
      if not isinstance(name, list) :
         name = [name]
      # populate the mask
      for w in what :
         for n in name :
            un.append(w, n)

      default = [map(lambda x:[], range(0,len(to)))]
      reqWeak = filter(lambda x: x in self.getWeakTables(), to)
      if any(map(lambda x: x in self.getWeakTables(), what)) :
         weak = True

      # if any asked table is weak but 'weak=False' then first get
      # masks without weak, then with weak and add the weak tables.

      # search the db
      try : 
         masks = self.getMask(un, weak=weak, wildcards=wildcards)
         # if no hits but we were not asked to search for weak
         # identifiers, try weak anyway
         if len(masks) == 0 and not external and not weak :
            masks = self.getMask(un, weak=True, wildcards=wildcards)
         if len(masks) > 0 and not weak and reqWeak:
            weak_masks = self.getMask(un, True, wildcards=wildcards)   
            for m in masks :
               for w in weak_masks:
                  if m.getIdentifiers("_id")[0] == w.getIdentifiers("_id")[0] :
                     for r in reqWeak:
                        if w.hasTable(r) and not m.hasTable(r):
                           m.copyTable(w, r)

      except NameError :
         return(default)
      # found tables
      found = map(lambda x:x.hasTable(to), masks)

      # generator for 
      def iden(tab, mzk) :
         for t in tab:
            yield(map(lambda ma:ma.getIdentifiers(t), mzk)[0])

      # make sure we do not query for something we do not know how to find
      if not any(map(lambda x:x in query.KNOWNTABLES, to)) :
         external = False

      # either found or explicitly told not to go online
      #if all(found) and found or not external or all([masks, to[0] == 'mask']):
      if all(found) and found or not external:
          if outmask :
             for ma in masks :
                for tab in ma.getTables():
                   if not tab in to:
                      ma.delTable(tab)
             return(masks)
          else :
             if masks :
                answer = []
                for ma in masks :
                   answer.append(list(iden(to, [ma])))
                return(answer)            
             else :
                return(default)
      
      # attach original search mask as well incase we had nothing
      # in the db
      masks.append(un)
      resmasks = []
      i = 1
      for ma in masks:
          if self.debug :
              print "#COMMENT No hit, trying external resources..",
          externalMask = query.fetch(self, ma, internal=len(masks) > 1, to=to)
          if not externalMask:
              if self.debug:
                  print "no hit, giving up"
              continue
          if self.debug :
              print "found it! "
          if learn:
              if self.debug:
                  print "#COMMENT adding the information"
              map(self.createIdTable, externalMask.getTables())
              self.setMask(externalMask)
          
          resmasks.append(externalMask)
          if i == len(masks) -  1 :
              break
          i = i + 1

      # return what we got
      if to[0] == 'mask':
          return(resmasks)
      else :
         answer = []
         for ma in resmasks :
            answer.append(list(iden(to, [ma])))
         return(answer)
      
   def getAllMmids(self) :
      """
      get all _ids
      """
      # get all _id
      self.c.execute('SELECT _id FROM _id')
      return(map(lambda x:x[0], self.c.fetchall()))

   def export (self, tables, weak=False) :
      """print all contents of the given tables to stdout

      Parameters:
      
      - `tables`: a list with the desired tables
      """
      
      mmids = self.getAllMmids()
      if tables[0] == "ALL" :
         tables = self.getIdTables()

      if any(map(lambda x: x in self.getWeakTables(), tables)) :
         weak = True
      
      # print header
      tmp = "\"" + "\",\"".join(tables) + "\""
      print tmp
      for mm in mmids :
         ma = mask({})
         ma.append('_id', mm)
         ma = self.getMask(ma, weak=weak)[0]
         if any(map(lambda x:x in tables, ma.getTables())) :
            i = 0
            # one mask one row
            for tab in tables :
               i = i + 1
               if ma.hasTable(tab) :
                  ide = ma.getIdentifiers(tab)
                  ide = map(str, ide)
                  tmp = "\"" + "\"|\"".join(ide) + "\""
                  sys.stdout.write(tmp)
               if i != len(tables) :
                  sys.stdout.write(",")
            # row is finished
            sys.stdout.write('\n')

   def stats (self, more=False) :
      """ print out basic statistics about the current database 
      """
      try:
         tabs = self.getIdTables()
         print "metmask db@" + self.db
         print "Known identifiers:"
         for t in tabs :
            self.c.execute("SELECT count(distinct " +t + ") FROM " + t)
            lines = map(lambda x:x[0], self.c.fetchall())
            print t + ":"
            print str(lines) + " rows"
         if more :
            self.c.execute("SELECT * FROM sources")
            print "Sources:"
            print '%-5s %-10s' % ('Code', 'Name')
            tmp = self.c.fetchall()
            for t in tmp :
                print '%-5s %-10s' % (str(t[0]),str(t[1]))
            self.c.execute("SELECT * FROM confidence")
            print "Confidence codes:"
            print '%-5s %-15s' % ('Code','Name')
            tmp = self.c.fetchall()
            for t in tmp :
                print '%-5s %-15s' % (str(t[0]),str(t[1]))
      except :
         raise dbError, "Could not perform a simple select statement"
         

         
