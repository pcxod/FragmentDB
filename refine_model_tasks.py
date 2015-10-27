'''
Created on 23.10.2015

@author: daniel
'''
import pprint



class Refmod(object):
    '''
    handles the refinement model of olex2
    '''

    def __init__(self):
      '''
      Constructor
      
      '''
      lstfile = 'd:\Programme\DSR\example\p21c.lst'
      final = False
      disag = False
      disargeelist = []
      with open(lstfile, 'r') as f:
        for num, line in enumerate(f):
          if line.startswith(" Final Structure Factor"):
            final = True
            strucfactline = num
          if line.startswith(' Total number of'):
            parameters = line.split()[6]
          if final and line.startswith(' wR2 ='):
            data = line.split()[7]
          if final and line.startswith(' Disagreeable restraints'):
            disag = True
          if disag: 
            disargeelist.append(line.split())
          if line.startswith(' Summary of restraints'):
            final = False
            disag = False         
      del disargeelist[:4]
      del disargeelist[-3:]
      print('residuals:', strucfactline, 'data:', data, 'parameters:', parameters)
      pprint.pprint(disargeelist)
    
    def open_lst_file(self):
      '''
      opens the resulting shelx lst file
      '''
      pass
      #ed = olx.GetVar('defeditor')
      #print(ed)
    
    def get_bond_list(self):
      pass
    
if __name__ == '__main__':
  ref = Refmod()