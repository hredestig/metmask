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

KNOWNTABLES = ['kegg', 'cas', 'cid', 'chebi', 'inchi', 'synonym']

MAX_SYN = 10
MAX_RET = 10


def nodecontents(nodes):
    """generator for getting the contents of an xml node 

    Parameters :
    -`nodes`: a minidom node
    """
    for n in nodes:
        for val in n.childNodes:
            yield val.nodeValue


def fetch(mm, un, internal, to):
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
    # if any([x in 'kegg' for x in to]):
    #     masks.extend(getKegg(mm, un))

    # if any([x in ['synonym', 'chebi', 'cid', 'cas', 'inchi'] for x in to]):
    #     masks.extend(getPubchem(mm, un))

    # TODO
    # this should somehow be done better so that obtained masks are
    # 'clustered' themselves before being checked with the database

    for pcmask in masks:
        suggestion = un.resolve(pcmask, mm)
        # check if we are allowed to merge and if, then do so
        # other wise don't
        if suggestion == 'merge' or not internal:
            return (pcmask)
    # means that we didn't merge
    return (False)
