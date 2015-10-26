'''
Created on 23.10.2015

@author: daniel
'''
from __future__ import print_function
import olex_core
from olexFunctions import OlexFunctions
import olx

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
      print('coord:---')
      print(olx.xf.au.GetAtomCrd('5'))
      print(olx.xf.au.GetCell())
      print('####################--------------###########')
      model = olex_core.GetRefinementModel()
      asym_unit = model['aunit']
      pf = pprint.pformat(model)
      #print(pf)
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
    def find_atoms_to_replace(self, frag_atoms, cell, remdist=1.2, only_this=None):
      '''
      this method looks around every atom of the fitted fragment and removes 
      atoms that are near a certain distance to improve the replace mode

      :param frag_atoms: atoms of the fitting fragment ['C1', '1', 'x', 'y', 'z']
      :param cell: unit cell parameters (list)
      :param remdist: distance below atoms shoud be deleted
      '''
      atoms_to_delete = []
      frag_coords = self.get_atomcoordinates(frag_atoms)
      atoms = self._residues
      for i in atoms:
        suffix = ''
        if i != 0:
          suffix = '_{}'.format(i)
        # i is the resideue number
        for y in atoms[i]:
          if y[0].startswith('Q'):
            # ignore q peaks:
            continue
          if only_this and only_this != y[5]:
            continue
          # y[4] is the part number
          if int(y[4]) == 0:
            for name in frag_coords:
              # name is the atom name
              if name == y[0]:
                # do not delete the fitted fragment
                continue
              at1 = [float(x) for x in frag_coords[name]]
              at2 = [float(x) for x in y[1]]
              resinum1 = self.get_atoms_resinumber(name)
              resinum2 = self.get_atoms_resinumber(y[0]+suffix)
              if at1 == at2 and resinum1 == resinum2:
                # do not delete atoms on exactly the same position
                # and same residue
                continue
              d = atomic_distance(at1, at2, cell)
              # now get the atom types of the pair atoms and with that
              # the covalence radius. 
              if d < remdist:
                atoms_to_delete.append(y[0]+suffix) 
      return sorted(atoms_to_delete)
