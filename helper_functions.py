'''
Created on 03.05.2015

@author: daniel

here are only functions that are completely independent of olex
'''
from FragmentDB_handler import FragmentTable

def call_profile(dbfile):
  import cProfile
  import pstats
  cp = cProfile.Profile()
  db=FragmentTable(dbfile)
  cp.runcall(db._get_fragment, 17)
  pstats.Stats(cp).strip_dirs().sort_stats('time').print_stats(20)
  

if __name__ == '__main__':
  dbfile = "fragment-database.sqlite"
  call_profile(dbfile)