'''
Created on 23.10.2015

@author: daniel
'''
from __future__ import print_function
import olex_core
from olexFunctions import OlexFunctions

import pprint

OV = OlexFunctions()

class Refmod(object):
    '''
    handles the refinement model of olex2
    '''

    def __init__(self):
      '''
      Constructor
      '''
      print('#######################refmodel######################')
        
    def get_bond_list(self):
      pass
    
    def get_atoms_list(self):
      model = olex_core.GetRefinementModel()
      asym_unit = model['aunit']
      pf = pprint.pformat(model)
      print(pf)
      print('#'*40)
      atoms = {}
      for residue in asym_unit['residues']:
        for atom in residue['atoms']:
          atoms[atom['aunit_id']] = [ atom['label'], atom['crd'][0], atom['part'] ]
      print(atoms)
      print(atoms[4])
      '''
      return
      for restr in model['sadi']:
        print(restr['value'], '#')
        for atom in restr['atoms']: # (atom_id, eqiv)
          print(atoms[atom[0]]['label'], '##')
          if atom[1]:
            print('$' + str(atom[1]), '###')
      '''

