import re

from metmask.mask import mask


class parser:
    def __init__(self, parent):
        parent.tables = ['kegg', 'synonym', 'cas', 'mpimp', 'mpimpclass']
        self.parent = parent

    def process(self):
        """ parse the Golm Metabolome Database standards library.
        such as
        those found at
        http://csbdb.mpimp-golm.mpg.de/csbdb/gmd/msri/gmd_msri.html
        """
        parent = self.parent
        parent.mm.setTableWeak('mpimpclass')
        # make sure that the necessary tables are in the db
        un = mask({}, parent.mm.idpatterns)
        l = parent.getLine()
        while l:
            # check and send off gathered data
            if re.match('^Name:', l):
                if un.nid(parent.mm) > 1:
                    parent.setMask(un)
                un = mask({}, parent.mm.idpatterns)
                x = re.findall('Name: (.+)', l)
                mpimpid = x[0].strip()
                un.append('mpimp', mpimpid, parent.confid,
                          parent.sourceid)
                # get synonyms from the title, the other synonyms are
                # just to full of crap
                spam = mpimpid.split('_')
                if spam[5].strip() != '':
                    un.append('synonym', spam[5], parent.confid,
                              parent.sourceid)
                    if re.match(".*\([0-9]+[tT][mM][sS]\)$",
                                spam[5]):
                        syn = re.findall("(.*)\([0-9]+[tT][mM][sS]\)$",
                                         spam[5])[0].strip()
                        un.append('synonym', syn,
                                  parent.confid, parent.sourceid)

            # ***** gather info *****
            # mpimp class 
            x = re.findall('Synonym: Metabolite \(class\):(.+)', l)
            if x:
                un.append('mpimpclass', x[0].strip(),
                          parent.mm.confidence['weak'],
                          parent.sourceid)
            # cas 1
            x = re.findall('CASNO: (.+)', l)
            if x:
                un.append('cas', x[0].strip(), parent.confid,
                          parent.sourceid)
            # cas 2
            x = re.findall('CAS\|(\d{1,7}-\d{2}-\d{1})_\(.+\)', l)
            if x:
                un.append('cas', x[0].strip(), parent.confid,
                          parent.sourceid)
            # kegg
            x = re.findall('KEGG\|([A-Z][0-9]{5})_.*\(.+\)', l)
            if x:
                un.append('kegg', x[0].strip(), parent.confid,
                          parent.sourceid)
            l = parent.getLine()
        # send of the last mask if we have one
        if un.nid(parent.mm) > 1:
            parent.setMask(un)
