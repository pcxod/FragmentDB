'''
Created on 09.10.2014

@author: Daniel Kratzert
'''
import sqlite3


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


class DatabaseRequest(object):
    def __init__(self):
        '''
        creates a database request
        '''
        self.con = sqlite3.connect(dbfile)

    def db_request(self, request):
        with self.con:
            self.cur = self.con.cursor()
            self.cur.execute(request)
            rows = self.cur.fetchall()
        return rows


class FragmentTable(DatabaseRequest, object):
    def __init__(self, dbfile):
        super(FragmentTable, self).__init__()

    def get_fragment_names(self):
            req = '''SELECT fragment.id, fragment.name FROM fragment ORDER BY name'''
            rows = self.db_request(req)
            for name in rows:
                print(name[0], name[1])


    def get_fragment(self, fragment_id):
        '''
        returns a full fragment with all atoms, restraints, ... as a dictionary
        :param fragment_id: id of the fragment in the database
        :type fragment_id: integer
        '''


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


class Restraints(object):
    def __init__(self):
        pass

    def get_restraints_from_fragment(self, fragment_name):
        pass

class Atom(object):
    def __init__(self):
        pass

    def get_fragment_atoms(self, fragment_name):
        pass


if __name__ == '__main__':
    dbfile = 'F:\GitHub\DSR-db\dk-database.sqlite'
    db = FragmentTable(dbfile)
    db.get_fragment_names()
    print('##################')
    found = db.find_fragment_by_name('CF3')
    print(found)

