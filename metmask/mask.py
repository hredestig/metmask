"""
How to use this module
======================

1. create a new mask with e.g. ``mask({})`` 
2. then append new identifiers with ``mask.append(tab, ide, ass, conf, src)``
3. use the provided getters and setters to set and retrieve identifiers from the mask.


Copyright (C) Henning Redestig
2009
See COPYING.txt for licensing details.
"""

import re
import sys

import metmask


class badMaskError(Exception):
    """Exception that should be raised when trying to form mask which is
    consisting of conflicting identfiers, eg H2O and O2
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class idMisMatchError(Exception):
    """Exception that should be raised when trying to insert append an
    identifier to table that does not match the legal pattern for that
    table.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def fixchebi(chebi):
    return str(chebi).replace('CHEBI:', '')


def fixcas(cas, constraints):
    """ try to fix a cas identifier that is on the integer form to its text
    form, i.e. 123456 becomes 123-45-6 and 123-45-6 is unchanged

    Parameters:
    - `cas` : string representing the cas number
    """
    inputcas = cas
    if 'cas' not in constraints or len(cas) < 4:
        raise idMisMatchError("Do not know how to fix : " + str(cas))
    if not re.match(constraints['cas'], cas):
        cas = cas[0:-3] + '-' + cas[-3:-1] + '-' + cas[-1]
    if not re.match(constraints['cas'], cas):
        raise idMisMatchError("not a CAS number: " + str(inputcas))
    return cas


def formula2dic(formula):
    """ take a chemical formula and turn it in to a dictionary with atom
    -> counts, ignore protons
    """
    vec = re.findall("[A-Z]{1}[a-z]{0,1}[0-9]*", formula)
    res = {}
    for v in vec:
        atom = re.findall("[A-Z]{1}[a-z]{0,1}", v)[0]
        numb = re.sub(atom, "", v)
        if not numb:
            numb = 1
            # ignore protons, no easy way of saying if they matter
            # only add the first atom, later may be complex waters
        if (not atom == 'H') & (atom not in res):
            res[atom] = int(numb)
    return res


def guessTable(string, constraints, mm=None):
    """ Return a vector with table names that match a set of given
    constraints. Only guesses tables that are already in the database
    if a pointer to the database object is given.

    Parameters:
    - `string`: a string, the identifiers at hand
    - `constraints` : a dictionary on the form tablename:regexp, if
      string matches regexp, the table name will be appended to the
      list to return.
    - `mm` : an instance of the metmask.dbi.db object
    """
    existingTables = set([])
    if mm:
        existingTables = set(mm.getIdTables())
    guesses = []
    string = str(string)
    for k in list(constraints.keys()):
        ok = True
        if existingTables:
            if k not in existingTables:
                ok = False
        if re.match(constraints[k], string) and ok:
            guesses.append(k)
    if not guesses:
        return [constraints['default']]
    return guesses


class mask:
    """ An object for managing a group of chemical identifiers that is
    assumed to be referring to a common parent chemical.
    """

    def __init__(self, cat=None, constraints=None):
        """ Initialize a `mask` object.

        Parameters:
        -`cat`: dictionary with the identifiers. It is formatted as
         {table:{identifier:[[associations],[source],[confidence]]}},
         where associations is the internal identifier of the
         association in the metmask db. The dictionary should normally
         not be created directly but by using the provided `append`
         function.
        `constraints`: dictionary with constraints specifying a
         regular expression for each table that the identifiers in
         that table must match.
        """
        self.cat = {} if cat is None else cat
        """ `dictionary` with the tables and identifiers that this mask is mapped to"""
        self.constraints = {} if constraints is None else constraints
        """ `dictionary` with the constraints"""
        self.caseSensitiveTables = ['inchi', 'formula', 'smiles', 'inchikey']
        """ tables that will not be down-cased
        """

    # low level functions
    def isEmpty(self):
        """ true if this maks does not map to anything """
        return (self.cat == {})

    def getTables(self):
        """ get a list with the tables specified for this mask """
        return (list(self.cat.keys()))

    def hasTable(self, table):
        """ true if the mask has all the specified tables.

         Parameters:
        - `table` : string or list specifiying the table(s)
        """
        if not isinstance(table, list):
            table = [table]
        res = [x in self.cat for x in table]
        return (all(res))

    def getAssoc(self, table, identifier):
        """ get the association identifiers for a certain table and identifier

        Parameters:
        - `table` : string specifiying the table
        - `identifier` : string specifiying the identifier
        """
        return (self.cat[table][identifier][0])

    def subsetToSource(self, src):
        """ get a new mask which only contains associations coming
        from a certain source
        Parameters:
        `src`:the wanted source
        """
        newmask = mask({})
        tables = self.getTables()
        for tab in tables:
            identifiers = self.getIdentifiers(tab)
            for ide in identifiers:
                ass = self.getAssoc(tab, ide)
                cnf = self.getConfidence(tab, ide)
                srz = self.getSource(tab, ide)
                for i in range(0, len(ass)):
                    if srz[i] == src:
                        newmask.append(tab, ide, cnf[i], srz[i], ass[i])
        return (newmask)

    def setAllAssoc(self, ass):
        """ set all associations in the associations in this mask to the
        provided association code
        Parameters:
        -`ass`:association code
        """
        tables = self.getTables()
        for tab in tables:
            identifiers = self.getIdentifiers(tab)
            for ide in identifiers:
                conf = self.getConfidence(tab, ide)
                self.cat[tab][ide][0] = [ass for x in range(0, len(conf))]

    def getConfidence(self, table, identifier):
        """ get the confidence identifiers for a certain table and identifier

        Parameters:
        - `table` : string specifiying the table
        - `identifier` : string specifiying the identifier
        """
        return (self.cat[table][identifier][1])

    def weakenIdentifiers(self, table, identifier):
        """ set the chosen identifier to be weak

        Parameters:
        - `table` : string specifiying the table
        - `identifier` : string specifiying the identifier
        """
        for i in range(0, len(self.cat[table][identifier][1])):
            self.cat[table][identifier][1][i] = metmask.WEAK_CONF

    def getSource(self, table=None, identifier=None):
        """ get the source identifiers for a certain table and identifier

        Parameters:
        - `table` : string specifiying the table, for all in unspecified
        - `identifier` : string specifiying the identifier,
           for all in unspecied
        """
        if not table and not identifier:
            src = []
            for tab in self.getTables():
                for ide in self.getIdentifiers(tab):
                    src.extend(self.getSource(tab, ide))
            return (list(set(src)))
        if not self.hasTable(table):
            return (None)
        if not self.hasId(table, identifier):
            return (None)
        return (self.cat[table][identifier][2])

    def getIdentifiers(self, table, weak=True):
        """ get the identifiers for a certain table

        Parameters:
        - `table` : string specifiying the table
        - `weak` : get weak identifiers as well
        """
        if not self.hasTable(table):
            return (None)
        ids = list(self.cat[table].keys())
        if not weak:
            ids = [x for x in ids if any([y != metmask.WEAK_CONF for y in self.getConfidence(table, x)])]
        return (ids)

    def copyTable(self, un, table):
        """ copy table from a mask to this mask

        Parameters:
        -`un` :a mask
        -`table` : the table to copy
        """
        self.cat[table] = un.cat[table]

    def delWeak(self, table):
        """ delete all weak identifiers in table
         Parameters:
        -`table` : desired table
        """
        ids = self.getIdentifiers(table, True)
        for id in ids:
            if all([x == metmask.WEAK_CONF for x in self.getConfidence(table, id)]):
                self.delIdentifiers(table, id)

    def delIdentifiers(self, table, identifier):
        """ drop the specified identifier

        Parameters:
        - `table` : string specifiying the table
        - `identifier` : string specifiying the identifier
        """
        del self.cat[table][identifier]
        if len(self.getIdentifiers(table)) == 0:
            self.delTable(table)

    def delTable(self, table):
        """ drop the specified table

        Parameters:
        - `table` : string specifiying the table
        """
        del self.cat[table]

    def hasId(self, table, identifier):
        """ true if the table has the given identifier in the specified table

        Parameters:
        - `table` : string specifiying the table
        - `identifier` : string specifiying the identifier
        """
        if not self.hasTable(table):
            return (False)
        return (identifier in self.cat[table])

    def toBIP(self, out, mm):
        """ obtain a graph representation of this mask
        -`source2master`: a dictionary defining which defines which the master
          identifier is for each source
        """
        fran = []
        till = []
        kall = []
        weak = []
        tables = self.getTables()
        masters = [mm.sourceid2master([src])[0] for src in self.getSource()]

        sources = [src for src in self.getSource() if mm.sourceid2master([src])[0] not in ['unknown', '_id']]
        if not sources:
            raise Exception("no sources have their master id set, re-import data and set master ids")
        for src in sources:
            rootIds = []
            srcmask = self.subsetToSource(src)
            if not mm.sourceid2master([src])[0] in tables:
                continue
            for ide in srcmask.getIdentifiers(mm.sourceid2master([src])[0]):
                if src in srcmask.getSource(mm.sourceid2master([src])[0], ide):
                    rootIds.append(ide)
            for rtId in rootIds:
                for tab in tables:
                    identifiers = []
                    if not srcmask.hasTable(tab):
                        continue
                    for ide in srcmask.getIdentifiers(tab):
                        if src in srcmask.getSource(tab, ide):
                            identifiers.append(ide)
                    for ide in identifiers:
                        if ide == rtId:
                            continue
                        if len([x for x in srcmask.getAssoc(str(mm.sourceid2master([src])[0]), rtId) if x in srcmask.getAssoc(tab, ide)]) > 0:
                            fran.append(str(mm.sourceid2master([src])[0]) + ":" + str(rtId))
                            till.append(str(tab) + ":" + str(ide).replace(",", ""))
                            kall.append(mm.sourceid2source([src])[0])
                            weak.append(mm.confidence['weak'] == \
                                        srcmask.getConfidence(tab, ide)[srcmask.getSource(tab, ide).index(src)])
        if len(fran) == 0:
            raise Exception("empty graph")
        for i in range(0, len(fran)):
            print("" + "\t".join(map(str, [str(fran[i]).replace('\t', ' '), \
                                         str(till[i]).replace('\t', ' '), \
                                         str(kall[i]).replace('\t', ' '), \
                                         str(weak[i]).replace('\t', ' ')])), file=out)

    def append(self, table, ident, cnf=None, src=None, ass=None):
        """ Add an identifier to this mask.

        Parameters:
        - `table`: string with the table type of this identifier
        - `ident`: string giving the identifier
        - `cnf` : the confidence code for this association
        - `src` : the source of evidence for this addition
        """

        # remove trailing whitespaces
        if isinstance(table, str):
            table = table.strip()
        if isinstance(ident, str):
            ident = ident.strip()
        if table == 'chebi':
            ident = fixchebi(ident)
        # check if we are adding a new sum formula -> not allowed
        if table == 'formula' and self.hasTable('formula'):
            curFormulas = self.getIdentifiers('formula')
            for cu in curFormulas:
                if formula2dic(cu) != formula2dic(ident):
                    print("skipping conspicuous mask: " + \
                                         str(self.getIdentifiers("_id")) \
                                         + " " + str(ident) + " and " + str(cu), file=sys.stderr)
                    return (None)
        # make identifiers always lower case and not a formula
        if 'encode' in dir(ident) and not table in self.caseSensitiveTables:
            ident = ident.lower()
        # chuck any pipes in the identifiers
        if re.match("\|", str(ident)):
            ident = re.sub("\|", "", str(ident))
        # check legality
        if table in self.constraints:
            if not table in guessTable(ident, self.constraints):
                # try to fix the identifier
                if table == 'cas':
                    ident = fixcas(ident, self.constraints)
                if not table in guessTable(ident, self.constraints):
                    raise idMisMatchError
                    # we are good, append:
        if not self.hasTable(table):
            self.cat[table] = {}
        if not self.hasId(table, ident):
            self.cat[table][ident] = [[], [], []]
        self.cat[table][ident][0].extend([ass])
        self.cat[table][ident][1].extend([cnf])
        self.cat[table][ident][2].extend([src])

    # high level functions, following functions must not touch 'cat' directly
    def nid(self, mm=None):
        """ return the number of tables in this mask
        Parameters:
        -`mm` : the dbi object, if given, strong tables are nout counted
        """
        tab = self.getTables()
        if mm:
            tab = [tab for tab in tab if tab not in mm.getWeakTables()]
        return (len(tab))

    def resolve(self, other, mm):
        """
        Returns 'prune' or 'merge' or None depending on if this mask
        is compatible with the other mask. Compatibility is checked by
        examining the following criteria (in this order):

        - If the two masks have no identifiers in common then return
          None.
         - If the two masks have non-equal formulas then return None
         - If both masks contain the nevermerge confidence level then
          return None.
         - If either mask contain the alwaysmerge confidence level then
          return 'merge'.

        - If the number of tables with overlapping identifiers are
          more than `MIN_OVERLAP`, and same confidence level then the
          return 'merge'

          + tables are only overlapping if they share non-weak identifiers

        - If a confidence level (excluding the nevermerge attribute)
          for a shared identifier is present in one mask but not the
          other, then return 'prune'
         - Otherwise, return 'merge'
         Parameters:
        - `other`: another instance of the `mask` class.
        - `mm`: an instance of the `db` class (compatibility is a
          database question).

        """
        tabSelf = self.getTables()
        tabOther = other.getTables()
        allTables = set(tabSelf + tabOther)
        overlap = []
        cself = []
        cother = []
        ideSelf = []
        ideOther = []
        diffConf = False
        if self.hasTable('formula') and other.hasTable('formula'):
            if formula2dic(self.getIdentifiers('formula')[0]) != formula2dic(other.getIdentifiers('formula')[0]):
                return (None)
        for tab in allTables:
            allIdentifiers = []
            if self.hasTable(tab):
                ideSelf = ideSelf + self.getIdentifiers(tab)
                allIdentifiers = allIdentifiers + self.getIdentifiers(tab)
            if other.hasTable(tab):
                ideOther = ideOther + other.getIdentifiers(tab)
                allIdentifiers = allIdentifiers + other.getIdentifiers(tab)
            allIdentifiers = set(allIdentifiers)
            for ide in allIdentifiers:
                if self.hasId(tab, ide):
                    cself = cself + self.getConfidence(tab, ide)
                if other.hasId(tab, ide):
                    cother = cother + other.getConfidence(tab, ide)
                if other.hasId(tab, ide) and self.hasId(tab, ide):
                    cS = self.getConfidence(tab, ide)
                    cO = other.getConfidence(tab, ide)
                    if (not all([x == mm.confidence['weak'] for x in cO]) and \
                                not all([x == mm.confidence['weak'] for x in cS])):
                        overlap = overlap + [tab]
                        # overlapping on different confidence levels -> False
                        diff = set(cS).difference(set(cO))
                        # potentially add more 'special' cofidence here
                        diff = diff.difference(set([mm.confidence['nevermerge']]))
                        if len(diff) > 0:
                            diffConf = True
                        #        if self.getIdentifiers('cas')[0] == '56-41-7' :
                        #          pdb.set_trace()
        numOverlap = len(set(overlap))

        # no overlap -> False
        if numOverlap == 0:
            return (None)
        # both tagged as nevermerge -> False
        elif mm.confidence['nevermerge'] in set(cother) and \
                        mm.confidence['nevermerge'] in set(cself):
            return (None)
        # fulfilled minimum overlap -> True
        elif numOverlap >= mm.MIN_OVERLAP:
            return ('merge')
        elif diffConf:
            return ('prune')
        # otherwise
        else:
            return ('merge')

    def merge(self, other):
        """ fuse a mask to this mask so that all association in the other
        masks are also in this mask.
         Parameters:

        - `other` : a differnet instance of `mask`
        """
        tables = other.getTables()
        for tab in tables:
            identifiers = other.getIdentifiers(tab)
            for ide in identifiers:
                conf = other.getConfidence(tab, ide)
                src = other.getSource(tab, ide)
                ass = other.getAssoc(tab, ide)
                for i in range(0, len(conf)):
                    self.append(tab, ide, \
                                conf[i], \
                                src[i], \
                                ass[i])

    def show(self, confidence=True, source=True, all=False, max=5, mm=None):
        """ Print some basic information about this mask.
        """
        if self.isEmpty():
            print("Empty mask")
            return (None)
        if self.hasTable('_id'):
            print("o-o-o-o:")
            print('mask:')
            print(str(self.getIdentifiers('_id')))
        for tab in self.getTables():
            if tab != '_id':
                print(tab + ":")
                n = 0
                identifiers = self.getIdentifiers(tab)
                for ide in identifiers:
                    if n < max or all:
                        print(ide, end=' ')
                        if confidence or source:
                            conf = self.getConfidence(tab, ide)
                            sour = self.getSource(tab, ide)
                            if mm:
                                sour = mm.sourceid2source(sour)
                            for i in range(0, len(conf)):
                                if source:
                                    print(sour[i], end=' ')
                                if confidence:
                                    print(conf[i], end=' ')
                        print()
                        n = n + 1
                    elif n == max and not all:
                        if len(identifiers) - max != 0:
                            print("(" + str(len(identifiers) - max) \
                                  + "more identifer(s)" + ")")
                        break

    def intersect(self, other):
        """ Return a new mask carrying the intersecting information between
        this mask and another mask.
         Parameters:
        - `other`: another instance of the `mask` class.
        """
        newmask = mask({})
        # the tables common to both masks
        tables = [x for x in other.getTables() if x in self.getTables()]
        for tab in tables:
            identifiers = [x for x in other.getIdentifiers(tab) if x in self.getIdentifiers(tab)]
            for ide in identifiers:
                conf = other.getConfidence(tab, ide)
                src = other.getSource(tab, ide)
                ass = other.getAssoc(tab, ide)
                for i in range(0, len(conf)):
                    newmask.append(tab, ide, \
                                   conf[i], \
                                   src[i], \
                                   ass[i])
                conf = self.getConfidence(tab, ide)
                src = self.getSource(tab, ide)
                ass = self.getAssoc(tab, ide)
                for i in range(0, len(conf)):
                    newmask.append(tab, ide, \
                                   conf[i], \
                                   src[i], \
                                   ass[i])
        return (newmask)

    def isSubset(self, other):
        """ true if all the information in this mask is also given
        in the other mask.
         Parameters:
        -`other`: another instance of the `mask` class.
        """
        raise idMisMatchError("Banned function")
        tables = other.getTables()
        for tab in tables:
            identifiers = other.getIdentifiers(tab)
            for ide in identifiers:
                if not self.hasId(tab, ide):
                    return (False)
        return (True)

    def anyNeverMerge(self):
        """ returns true if this mask has an identifier marked as  nevermerge
        """
        tables = self.getTables()
        for tab in tables:
            identifiers = self.getIdentifiers(tab)
            for ide in identifiers:
                if any([x == metmask.NEVERMERGE_CONF for x in self.getConfidence(tab, ide)]):
                    return (True)
        return (False)

    def subtract(self, other, onlyTables=[]):
        """ subtract the content of the other mask from this mask.
         Parameters:
        -`other` : another instance of the `mask` class.
        -`onlyTables` : a list of tables that should be
         considered. Defaults to all tables if not specified.
        """
        # the tables common to both masks
        tables = [x for x in other.getTables() if x in self.getTables()]
        if onlyTables:
            tables = [x for x in tables if x in onlyTables]
        for tab in tables:
            identifiers = [x for x in other.getIdentifiers(tab) if x in self.getIdentifiers(tab)]
            for ide in identifiers:
                self.weakenIdentifiers(tab, ide)
                if tab == 'formula':
                    self.delIdentifiers(tab, ide)
