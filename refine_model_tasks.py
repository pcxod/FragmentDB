'''
Created on 23.10.2015

@author: daniel
'''
from __future__ import print_function
import olex_core
from olexFunctions import OlexFunctions
import olx

OV = OlexFunctions()

class Refmod(object):
    '''
    handles the refinement model of olex2
    '''

    def __init__(self):
      '''
      Constructor
      '''
    
    def open_lst_file(self):
      '''
      opens the resulting shelx lst file
      '''
      ed = olx.GetVar('defeditor')
      print(ed)
    
    def get_bond_list(self):
      pass
    
    