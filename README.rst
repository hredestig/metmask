Installation
============

Old legacy project, kept here for reference purpose only.

This is a python package for managing lists of chemical identifiers
with focus on metabolomics. You need to have the python interpreter
installed in order to use it. Python can be obtained at
http://www.python.org. You also need 'SOAPpy' installed (which depends
on 'fpconst').

Metmask is a command line program and therefore does not offer any
normal `point-and-click' interface. Instead you need to use the
command line via e.g. a graphical shell like the Gnome Terminal. 

Metmask is currently developed for Python 2.5 but also tested on 2.6

Linux (Debian/Ubuntu)
---------------------

Unpack::

  tar -zxvf metmask-X.Y.Z.tar.gz
  cd metmask-X.Y.Z

then either install for all users::

  sudo python setup.py install

or locally in your home folder::

  python setup.py install --home=~

Simple usage
============

Installing the metmask package should also put the ``metmask`` script
somewhere in your executable path. You can then run ``metmask`` simply
by calling it from your shell, e.g.::

  metmask --version

Metmask comes with an example database and you can get obtain basic
  statistics for by calling :: 

  metmask -s

You can query it by one identifier at a time, or by passing a file
containing your identifiers. Example for querying for a single ID
fetching KEGG and synonyms:
   
  metmask -q 77-92-9 -g kegg,synonym

You can also make your own database by first importing a list of
identifiers that may like this::

  synonym,cas,kegg
  water,7732-18-5,c00001
  alanine,2899-44-7|302-72-7,c00041|c00133

first line names the type of identifier comma separated. Following
lines contains your identifiers with multiple identifiers separated by
| (pipe) (though this can be configured). This example file should be
the metmask-X.Y.Z folder that you just unpacked. Import this to a new 
database by issuing::

  metmask -i test-identifiers.txt -p simple --db newdb

Where ``-p simple`` indicates that the ``simple`` parser should be
used to parse the file.

The created database can be augmented with data from e.g. ChEBI by
doing (this will take a while so please be patient) ::

  metmask -i chebi --db newdb 

When executing ``metmask`` without any switches the program will wait
for identifiers to query the database for from standard input. If you
want to query for a single identifier use the ``-a`` switch, e.g::

  metmask -q alanine -g formula,chebi --db newdb

If you have a file of identifiers that you want to query the database
for you can simply do e.g.::

  metmask -g chebi,kegg -d newdb < test-query.txt

and use standard command-line redirection to save the result. E.g::

  metmask -g chebi,kegg -d newdb < test-query.txt > result.txt

See::

  metmask -h

for a description of the available commands and the user manual
(metmask-X.Y.Z/doc/manual.txt) for further documentation.

