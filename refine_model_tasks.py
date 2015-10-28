'''
Created on 23.10.2015

@author: daniel
'''
import pprint
try:
  from olexFunctions import olx
  OV = OlexFunctions()
  IT = ImageTools()
except:
  pass



class Refmod(object):
    '''
    handles the refinement model of olex2
    '''

    def __init__(self):
      '''
      Constructor
      
      '''
      self.lstfile = 'd:\Programme\DSR\example\p21c-test.lst'
      #self.lstfile = '/Users/daniel/Documents/DSR/example/p21c.lst'
      
    def fileparser(self):
      '''
      gathers the residuals of the lst file
      '''  
      final = False
      disag = False
      disargeelist = []
      with open(self.lstfile, 'r') as f:
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
        return disargeelist
    

    
    def get_bond_list(self):
      pass
    
    
    def results_window(self):
      '''
      display results in a window
      '''
      pop_name = "Residuals"
      screen_height = int(olx.GetWindowSize('gl').split(',')[3])
      screen_width = int(olx.GetWindowSize('gl').split(',')[2])
      box_x = int(screen_width*0.1)
      box_y = int(screen_height*0.1)
      width, height = 500, 520
      html = r"""
            <table border="0" cellpadding="0" cellspacing="5" width="100%" > 
              {}
            </table
              """.format(targs)
      for line in disagreelist:
        td = r"""<td> {} </td>""".format(line)
        for r in td:
          row = r"""<tr> {}  </tr>""".format(r)
          
        
      OV.write_to_olex('large_fdb_image.htm', html)
      olx.Popup(pop_name, "large_fdb_image.htm",  b="tcrp", t="View Fragment", w=width,
                h=height, x=box_x, y=box_y)
    
if __name__ == '__main__':
  ref = Refmod()
  ref.fileparser()
  
  