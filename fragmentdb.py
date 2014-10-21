'''
Created on 09.10.2014

@author: Daniel Kratzert
'''

__metaclass__ = type  # use new-style classes
import sqlite3
import sys
from sqlite3 import OperationalError
print(sys.version)
print

def dice_coefficient(a, b):
    '''
    dice coefficient 2nt/na + nb
    Compares the similarity of a and b
    :param a: string
    :param b: string
    '''
    a = a.lower()
    b = b.lower()
    if not len(a) or not len(b): return 0.0
    if len(a) == 1:  a=a+'.'
    if len(b) == 1:  b=b+'.'
    a_bigram_list=[]
    for i in range(len(a)-1):
        a_bigram_list.append(a[i:i+2])
    b_bigram_list=[]
    for i in range(len(b)-1):
        b_bigram_list.append(b[i:i+2])
    a_bigrams = set(a_bigram_list)
    b_bigrams = set(b_bigram_list)
    overlap = len(a_bigrams & b_bigrams)
    dice_coeff = overlap * 2.0/(len(a_bigrams) + len(b_bigrams))
    dice_coeff = 1-dice_coeff # invert the result
    if dice_coeff < 0.5:  # make a cutoff for the best matches
        return 0.0
    return round(dice_coeff, 6)


class DatabaseRequest():
    def __init__(self, dbfile):
        '''
        creates a database request
        '''
        # open the database
        self.con = sqlite3.connect(dbfile)
        self.con.text_factory = str

    def db_request(self, request, *args):
        '''
        Performs a SQLite3 database request with "request" and optional arguments
        to insert parameters via "?" into the database request.
        :param request:
        :type request:
        '''
        with self.con:
            # set the database cursor
            self.cur = self.con.cursor()
            # make the request
            try:
                self.cur.execute(request, args)
            except OperationalError as e:
                print(e)
                return False
            # fetch all rows
            rows = self.cur.fetchall()

        return rows


class FragmentTable(DatabaseRequest):
    def __init__(self, dbfile):
        super(FragmentTable, self).__init__(dbfile)
        self.all_fragments = self.get_all_fragment_names()
        if not self.all_fragments:
            print('no database items found')
            sys.exit()
        self.all_names = [i[1].lower() for i in self.all_fragments]

    def __contains__(self, name):
        if type(name) == str():
            found = name.lower() in self.all_names
            return found

    def __iter__(self):
        return (i[1] for i in self.all_fragments)

    def get_all_fragment_names(self):
        '''
        returns all fragment names in the database, sorted by name
        '''
        req = '''SELECT fragment.id, fragment.name FROM fragment ORDER BY name'''
        rows = self.db_request(req)
        return rows

    
    def get_fragment(self, fragment_id):
        '''
        returns a full fragment with all atoms, restraints, atom types as a dictionary
        :param fragment_id: id of the fragment in the database
        :type fragment_id: integer
        '''
        req_atoms = '''SELECT atoms.name, atoms.element, atoms.x, atoms.y, atoms.z
            FROM fragment, atoms on fragment.id=atoms.fragmentid WHERE
            fragment.id = {}'''.format(fragment_id)
        req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
            FROM Restraints WHERE FragmentId = ?'''
        atomrows = self.db_request(req_atoms)
        restraintrows = self.db_request(req_restr, fragment_id)
        return (restraintrows, atomrows)


    def find_fragment_by_name(self, name, selection=5):
        '''
        find a fragment by its name in the database. This method will output a
        selection of (default=5) best hits.
        :param name: (part of) the name of a fragment to find
        :type name: string
        :return fragment_id: list of id numbers of the found fragments e.g. [1, 5]
        :type fragment_id: integer
        '''
        req = '''SELECT fragment.id, fragment.name FROM fragment'''
        frags = self.db_request(req)
        fragment_ids = self._search_name(name, frags, selection)
        return fragment_ids


    def _search_name(self, search_string, frags, selection=5):
        '''
        searches the names in the database for a given name
        :param search_string: search for this string
        :type search_string: string
        :param frags: fragments in the database
        :type frags: list
        :param selection: return this amount of results
        :type selection: integer
        '''
        search_results = {}
        for i in frags:
            db_entry = i[1]
            coefficient = dice_coefficient(search_string, db_entry)
            search_results[coefficient] = i
        # select the best [selection] results:
        selected_results = [search_results[i] for i in sorted(search_results)[0:selection]]
        return selected_results


    def store_fragment(self, fragment_name, atoms, restraints, tag, reference, comment):
        '''
        Store a new fragment into the database
        :param fragment_name: full chemical name of the fragment
        :type fragment_name: string
        :param atoms: atoms of the fragment ['name', 'atomic number', 'x', 'y', 'z']
        :type atoms: list
        :param restraints: [['DFIX', '1.564', 'C1 C2'], ['SADI', 'C2 C3 C4 C5']]
        :type restraints: list of list
        :param tag: short name tag (not mandatory)
        :type tag: string
        :param reference: source of the fragment
        :type reference: string
        :param comment: any comment
        :type comment: string
        '''
        pass


class Restraints(DatabaseRequest):
    def __init__(self, dbfile):
        super(Restraints, self).__init__(dbfile)


    def get_restraints_from_fragmentId(self, fragment_id):
        req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
            FROM Restraints WHERE FragmentId = ?'''
        restraintrows = self.db_request(req_restr, fragment_id)
        return(restraintrows)



if __name__ == '__main__':
    dbfile = 'F:\GitHub\DSR-db\dk-database.sqlite'
    db = FragmentTable(dbfile)
    print(db.get_all_fragment_names())
    #import cProfile
    #cProfile.run("allnames()", "foo.profile")
 #   res = Restraints(dbfile)
  #  for r in res.get_restraints_from_fragmentId(15):
  #      pass
        #print(r)

 #   print(hasattr(db, '__iter__'))
 #   if 'OC(CF3)3' in db:
 #       print('yes')
 #   else:
 #       print('no benz')

    # for i in db:
    #     print(i)

    #print('##################')
    #found = db.find_fragment_by_name('CF3', selection=3)
    #print(found)
    #for i in db.get_fragment(fragment_id=2):
    #    print(i)
