"""
This is the main interface to the metmask database. Metmask is a
tool for mapping chemical compounds as from biological samples as
identified by high-throughput chromatographic methods. Different
analytes are attempted to be mapped to the same original metabolite
so that e.g. two TMS derivates of the same metabolite can easily be
mapped to that original metabolite. Rare metabolites, such as
D-Amino acids or L-sugars should map to the more common L-Amino
acid and D-sugar forms as chromatographics methods rarely can
detect chiral isomers and hence are assumed to be just
mis-identified. Building is designed to require as little manual
intervention as possible.

Copyright (C) Henning Redestig
2009
See COPYING.txt for licensing details.
"""
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from builtins import input
from builtins import str
from builtins import range
from optparse import OptionParser
import fileinput, sys, urllib.request, urllib.error, urllib.parse, configparser, os, re, traceback
import metmask
import metmask.mask
import metmask.dbi
import metmask.parse
import socket
from metmask.parse import main as pmain


def oUrl(url):
    req = urllib.request.Request(url)
    try:
        handle = urllib.request.urlopen(req)
        next(handle)
        handle = urllib.request.urlopen(req)
    except:
        queryError("bad url, correct the metmask config file")
    return (handle)


def queryMessage(message):
    print()
    "#COMMENT:" + message


def queryWarning(message):
    print()
    "#COMMENT:" + message
    print("#COMMENT:", message, file=open(errlog, "a"))


def queryError(message, inst=None):
    """generic query error, used to die gracefully without printing ugly tracebacks
    """
    message = message + " ** python error: " + str(inst)
    print("#ERROR:" + message, file=sys.stderr)
    if inst and options.verbose:
        traceback.print_exc(file=open(errlog, "a"))
        raise inst
    else:
        try:
            quit()
        except:
            raise SystemExit


def doQuery(iden):
    """ query the database for an identifier. this function needs a major cleanup
    """
    nastring = "NA"
    if not options.noquote:
        nastring = "'NA'"
    inputIden = iden.strip()
    # if we get nothing then return just an empty line, this makes
    # scripting with empty lines easier
    if inputIden == "":
        print()
        return (0)

    try:
        iden = pmain.fixLine(iden, sep='\|')
    except ValueError:
        iden = iden.lower().strip().split(simple['sep2'])
    to = pmain.fixLine(options.goal, sep=simple['sep1'])

    if options.verbose:
        queryMessage("Querying metmask database..")
    if not options.table:
        table = []
        for ide in iden:
            table.extend(metmask.mask.guessTable(ide, mm.idpatterns, mm))
    else:
        table = options.table.split(simple['sep1'])
        # caution: tables should be interpreted as per identifier
        # regardless of length?
        #
        # if not len(table) == len(iden) : queryError("specify none
        # or exactly all input tables")
    if options.verbose:
        queryMessage("Looking for " + str(table) + " = " + str(iden))
    try:
        if iden[0] == 'NA':
            res = []
        else:
            res = mm.simpleQuery(name=iden, what=table, to=to,
                                 external=options.online,
                                 learn=options.universal,
                                 weak=options.universal,
                                 outmask=options.output in ['mask', 'graph'],
                                 wildcards=options.wild)
    except Exception as inst:
        queryError("unknown table maybe?", inst)
    if options.output == 'flat':
        output_flat(res, nastring, to)
    elif options.output == 'mask':
        output_mask(res)
    elif options.output == 'graph':
        output_graph(res)


def output_graph(res):
    for r in res:
        try:
            r.toBIP(out=sys.stdout, mm=mm)
        except Exception as inst:
            queryError("Graph to generate is empty, maybe not enough identifier types?", inst)


def output_mask(res):
    for r in res:
        r.show(all=True, mm=mm)


def output_flat(res, nastring, to):
    if not res:
        sys.stdout.write("\"" + "\",\"".join([nastring] * len(to)) + "\"\n")
    else:
        ran = list(range(0, len(res)))
        if options.onehit:
            ran = [0]
        for k in ran:
            subres = res[k]
            if not subres:
                sys.stdout.write(nastring)
            for i in range(0, len(subres)):
                if not options.noquote:
                    sys.stdout.write("'")
                if not subres[i]:
                    if not options.noquote:
                        sys.stdout.write("NA'")
                    else:
                        sys.stdout.write("NA")
                else:
                    ranID = list(range(0, len(subres[i])))
                    if options.first:
                        ranID = [0]
                    for j in ranID:
                        if not options.noquote:
                            sys.stdout.write('"%(hit)s"' % {'hit': subres[i][j]})
                        else:
                            sys.stdout.write('%(hit)s' % {'hit': subres[i][j]})
                        if j != max(ranID):
                            sys.stdout.write(simple['sep2'])
                        else:
                            if not options.noquote:
                                sys.stdout.write("'")
                if i != len(subres) - 1:
                    sys.stdout.write(simple['sep1'])
            sys.stdout.write("\n")


def quit():
    try:
        mm.close()
    except:
        pass
    raise SystemExit


if __name__ == '__main__':

    # windows safe hopefully
    try:
        if os.name == 'posix':
            configfile = os.path.join(os.getenv('HOME'), '.metmask.cfg')
            errlog = os.path.join(os.getenv('HOME'), 'metmask-error.log')
        else:
            configfile = metmask.dbi.determine_path() + '/.metmask.cfg'
            errlog = metmask.dbi.determine_path() + 'metmask-error.log'
    except:
        configfile = '.metmask.cfg'
        errlog = 'metmask-error.log'

    dbfile = metmask.dbi.determine_path() + "/data/metmask-db"

    # try to fetch the configurations, if they are not there, write the
    # default config file
    config = configparser.RawConfigParser()
    configurations = config.read(configfile)
    if not configurations:
        config = configparser.RawConfigParser()
        config.add_section('general')
        config.set('general', 'db', dbfile)
        config.add_section('simple')
        config.set('simple', 'na', '(^[nN][aA]$)|(-*$)')
        config.set('simple', 'sep1', ',')
        config.set('simple', 'sep2', '|')
        config.set('general', 'goal', 'synonym')
        config.set('general', 'token', '')
        config.set('general', 'confidence', 'good')
        config.set('general', 'ask', '0')
        config.set('general', 'minoverlap', '2')
        config.set('general', 'kegg', 'ftp://ftp.genome.jp/pub/kegg/ligand/compound/compound')
        config.set('general', 'cyc', 'ftp://ftp.plantcyc.org/Pathways/plantcyc_compounds.20091014')
        config.write(open(configfile, 'wb'))
        configurations = config.read(configfile)
    defaults = {}
    simple = {}
    # if there was a config file, fix the defaults
    if configurations:
        for pair in config.items('general'):
            defaults[pair[0]] = pair[1]
        for pair in config.items('simple'):
            simple[pair[0]] = pair[1]

    # if not defaults.has_key('token') :
    #    defaults['token'] = ''

    # set up the command line options parser
    parser = OptionParser(usage="""%prog [X] [options]
    X represents a stream of identifiers from STDIN
    or a single identifier if -a is set or a file with
    identifiers, one per line. If X is not provided one
    or more of the options below must be set.""",
                          version="This is " + metmask.__name__ + " v" + metmask.__version__)
    parser.add_option("-a", "-q", "--query", "--as-is-id", action="store", dest="asIs", type="string",
                      metavar="IDENTIFIER",
                      help="A single  identifier to query the data base for to be interpreted as is. Also searches "
                           "for weak (ambiguous) identifiers is the universal argument is set ")
    parser.add_option("-w", "--wild", action="store_true",
                      help="Input query may contain wild cards ('_' matches any single character and '%' matches any "
                           "sequence of characters). ",
                      dest="wild")
    parser.add_option("-t", "--table", action="store", type="string",
                      help="Identifier is of type TABLE (i.e. identifier type). Metmask tries to guess which table is "
                           "meant in case unspecified. Multiple tables can be given comma separated.",
                      metavar="TABLE",
                      dest="table")
    parser.add_option("-g", "--goal", action="store", type="string",
                      help="Fetch goal identifier of TABLE or use 'ALL' to get all entries. Multiple tables can be "
                           "given pipe separated.",
                      metavar="TABLE",
                      default=defaults['goal'],
                      dest="goal")
    parser.add_option("-d", "--db", action="store", type="string",
                      dest="db",
                      help="Use database located at PATH",
                      metavar="PATH",
                      default=defaults['db'])
    parser.add_option("-s", "--stats", action="store_true",
                      help="Print some statistics of the current database. Print more stats if the universal "
                           "argument is set.",
                      dest="stats")
    parser.add_option("-v", "--verbose", "--debug", action="store_true",
                      help="Be more verbose, print debug information",
                      dest="verbose")
    parser.add_option("-i", "--import", action="store", type="string",
                      help="Populate database using this file. For importing the KEGG compounds file or the PlantCyc "
                           "compounds file, FILE can also be one of the keywords 'cyc' or 'kegg' in order to use the "
                           "ftp provided files. Local files have precedence.\nUniversal argument: Never ask before "
                           "merging on input, guess what to do",
                      metavar="FILE",
                      dest="popFile")
    parser.add_option("-n", "--name", action="store", type="string",
                      help="Each import is given a name, usually this equals the filename or name of external "
                           "database but can be set explicitly with this option",
                      metavar="NAME",
                      dest="source")
    parser.add_option("--minoverlap", action="store", type="string",
                      default=defaults['minoverlap'],
                      help="The minimum number of tables in which two masks must overlap to be merged upon this import",
                      metavar="MINOVERLAP",
                      dest="minoverlap")
    parser.add_option("-c", "--confidence-code", action="store", type="string",
                      default=defaults['confidence'],
                      help="The confidence code of the data to import",
                      metavar="CODE",
                      dest="confidence")
    parser.add_option("-p", "--parser", action="store",
                      help="Use this parser to populate the database",
                      choices=[re.sub("^_", "", x) for x in metmask.parse.PARSERS],
                      dest="parser")
    parser.add_option("-o", "--output", action="store",
                      help="Output modes, flat: comma delimited output, mask: a text represenation of all info in the "
                           "mask, graph: a connection graph to use for visualization in e.g. cytoscape.",
                      default='flat', choices=['flat', 'mask', 'graph'],
                      dest="output")
    parser.add_option("-M", "--master", action="store", help="Master identifier for this import", dest="master")
    parser.add_option("-m", "--merge-interactively", action="store_true",
                      help="If set, ask about merging masks upon insert.",
                      dest="merge")
    parser.add_option("-r", "--remove", action="store", metavar="MMID",
                      help="Drop this identifier. Must be combined with the tabletype (-t switch). Can not be undone! "
                           "\n Universal argument: Do not ask about deletions before",
                      dest="delMask")
    parser.add_option("-x", "--export", action="store", metavar="TABLES",
                      help="Export all information in the provided tables (comma separated, potentially quoted). Use "
                           "the keyword ALL for exporting all known tables.",
                      dest="export")
    parser.add_option("-u", "-F", "--universal", action="store_true",
                      help="Universal argument. If set, causes some arguments to act differently.",
                      dest="universal")
    parser.add_option("-f", "--first", action="store_true",
                      help="Only return one identifier of each queried table. Suppress printing multiple identfiers.",
                      dest="first")
    parser.add_option("-e", "--external", action="store_true",
                      help="Try to query PubChem and KEGG if no suitable match could be found. Universal argument: "
                           "Save retrieved information in the local database",
                      dest="online")
    parser.add_option("-1", "--one-hit", action="store_true",
                      help="Only return the first hit  from the database",
                      dest="onehit")
    parser.add_option("-Q", "--no-quote", action="store_true",
                      help="Suppress addition of quotes to output",
                      dest="noquote")
    parser.add_option("-S", "--synchronize", action="store_true",
                      help="Upon import, minimize the creation of new masks, instead just populate the existing ones "
                           "and ignore all other input",
                      dest="boost")

    (options, args) = parser.parse_args()

    if options.merge:
        options.merge = defaults['ask']

    # get a connection to the db
    # not a new database but no write access
    try:
        mm = metmask.dbi.db(db=options.db, ask=options.merge, debug=options.verbose, minoverlap=options.minoverlap)
    except Exception as inst:
        queryError('bad db @ ' + options.db, inst)

    # get some statistics
    if options.stats:
        try:
            mm.stats(options.universal)
        except Exception as inst:
            queryError("", inst)
        quit()

    # delete a mask
    if options.delMask:
        if not options.table:
            queryError("specify identifier type (using the -t switch) when dropping identifiers")
        if options.universal:
            choice = 'yes'
        else:
            choice = input("really try to drop identifier " + str(options.delMask) + " (yes/no): ")
        while (choice not in ['yes', 'no']):
            choice = input("please answer yes or no: ")
        if choice == 'yes':
            un = metmask.mask.mask({})
            un.append(options.table, options.delMask)
            try:
                if options.table == "_id":
                    mm.dropMask(un)
                else:
                    mm.dropIdentifiers(un)
            except Exception as inst:
                queryError("maybe no such identifier", inst)
        quit()

    # export
    if options.export:
        tables = pmain.fixLine(options.export, sep=",")
        badTables = [x for x in tables if x not in mm.getIdTables() + ['ALL']]
        if badTables:
            queryError("unknown table" + str(badTables) + \
                       " choose one or more (comma-delimited) from " + \
                       str(mm.getIdTables()) + \
                       " or the keyword ALL for exporing all tables")
        mm.export(tables, weak=options.universal)
        quit()

    # populate the database
    elif options.popFile or options.parser:
        if options.popFile:
            if not options.source:
                options.source = os.path.basename(options.popFile)
                # cyc via ftp
        if options.popFile == 'cyc' and not os.path.isfile(options.popFile):
            options.popFile = oUrl(defaults['cyc'])
            options.parser = 'cyc'
            options.source = 'plantcyc.org'
            # kegg via ftp
        if options.popFile == 'kegg' and not os.path.isfile(options.popFile):
            socket.setdefaulttimeout(60)
            options.popFile = oUrl(defaults['kegg'])
            if not options.parser:
                options.parser = 'kegg'
                options.source = 'genome.jp'
                # chebi via http post
        if options.popFile == 'chebi':
            # popfile True means that the parser should use the db instead
            options.popFile = True
            if options.boost:
                queryWarning("synchronization option not needed, parser chebi always augment to current database")
                options.boost = False
            if not options.parser:
                options.parser = 'chebi'
                options.source = 'www.ebi.ac.uk/chebi'

        if options.popFile == 'chemspider' or options.parser == 'chemspider':
            if options.popFile == 'chemspider':
                options.popFile = True
            if options.boost:
                queryWarning("synchronization option not needed, parser chemspider always augment to current database")
                options.boost = False
            if not options.parser:
                options.parser = 'chemspider'
                options.source = 'chemspider.com'

        if options.popFile == 'pubchem' or options.parser == 'pubchem':
            if options.popFile == 'pubchem':
                options.popFile = True
            if options.boost:
                queryWarning("synchronization option not needed, parser pubchem always augment to current database")
                options.boost = False
            if not options.parser:
                options.parser = 'pubchem'
                options.source = 'pubchem.ncbi.nlm.nih.gov'

                # ..user provided a file or we die
        if not options.popFile and options.parser:
            queryError("Need both file and parser to populate the database")
        if options.popFile and options.parser and not options.source:
            queryError("Provide a name for this import using the -n switch")
        try:
            # parser = metmask.parse.getParser(options.parser)
            myimporter = pmain.importer(mm, parsertype=options.parser,
                                        fileobj=options.popFile,
                                        source=options.source,
                                        confidence=options.confidence,
                                        sep1=simple['sep1'],
                                        sep2=simple['sep2'],
                                        na=simple['na'],
                                        resolve=not options.universal,
                                        boost=options.boost,
                                        master=options.master,
                                        token=defaults['token'])
            myimporter.parser.process()
        except NameError as inst:
            queryError("Choose one of following parsers: " + str(metmask.parse.PARSERS), inst)
        except Exception as inst:
            queryError("", inst)
        quit()

    # single query
    elif options.asIs:
        doQuery(options.asIs)
        quit()

    # nothing else now query STDIN instead
    try:
        for iden in fileinput.input(args):
            doQuery(iden)
    except IOError as inst:
        queryError("Invalid choice of options or bad input", inst)

    quit()
