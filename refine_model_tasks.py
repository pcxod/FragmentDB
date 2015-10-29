'''
Created on 23.10.2015

@author: daniel
'''
#import pprint
from helper_functions import REL_RESTR_CARDS
try:
  import olx  # @UnresolvedImport
  from olexFunctions import OlexFunctions
  #from ImageTools import ImageTools
except:
  pass
try:
  OV = OlexFunctions()
  #IT = ImageTools()
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
      self.htm = html_Table()
      
    def fileparser(self, lstfile):
      '''
      gathers the residuals of the lst file
      '''  
      final = False
      disag = False
      disargeelist = []
      with open(lstfile, 'r') as f:
        for line in f:
          if not line.split():
            continue
          if line.startswith(" Final Structure Factor"):
            final = True
#          if line.startswith(' Total number of'):
#            parameters = line.split()[6]
#          if final and line.startswith(' wR2 ='):
#            data = line.split()[7]
            #disargeelist.append(['Results', 'data:', data, 'parameters:', parameters])
          if final and line.startswith(' Disagreeable restraints'):
            disag = True
            continue
          if line.startswith(' Summary of restraints'):
            final = False
            disag = False
          if disag: 
            fline = self.lineformatter(line.split())
            # fline is a list of list
            disargeelist.append(fline)
        return disargeelist
    
    def lineformatter(self, line):
      '''
      takes care of some extra things with different restraints. For example
      RIGU xy should be in one column
      '''
      if line[0].startswith('Observed'):
        return line
      for num, i in enumerate(line):
        if i[0].isalpha():
          pos = num
          break
      line[pos:] = [' '.join(line[pos:])]
      tline = ' '.join(line)
      for n in REL_RESTR_CARDS:
        if n in tline:
          line = ['-', '-']+line
          break
      return line
        
    def open_listfile(self, title = "Select a .lst file", 
                            ffilter = '*.lst; *.LST',
                            location = '',
                            default_name = '' ):
      lstfile = olx.FileOpen(title, ffilter, location, default_name)
      return lstfile 
    
    def results_window(self):
      '''
      display results in a window
      '''
      lstfile = self.open_listfile()
      if not lstfile:
        print('No file selected')
        return
      pop_name = "Residuals"
      screen_height = int(olx.GetWindowSize('gl').split(',')[3])
      screen_width = int(olx.GetWindowSize('gl').split(',')[2])
      box_x = int(screen_width*0.1)
      box_y = int(screen_height*0.1)
      width, height = 600, 520
      filedata = self.fileparser(lstfile)
      if not filedata:
        filedata =['']
      header = ['<h2>List of most disagreeable restraints</h2>']
      footer = ['<br><br> Use "MORE 4" to get an extensive list of all restraints. ']
      html = self.htm.table_maker(header, filedata, footer)
      OV.write_to_olex('large_fdb_image.htm', html)
      olx.Popup(pop_name, "large_fdb_image.htm",  b="tcrp", t="View Fragment", w=width,
                h=height, x=box_x, y=box_y)
    

class html_Table(object):
  '''
  html table generator
  '''
  def __init__(self):
    pass

  
  def table_maker(self, header=[''], tabledata=[''], footer=['']):
    '''
    builds a html table out of a datalist from the final 
    cycle summary of a shelxl list file.
    '''
    table=[]
    data = []
    for line in tabledata:
      table.append(self.row(line))
    data.extend(header)
    data.extend(table)
    data.extend(footer)
    html = r"""
      {0}
    <table border="0" cellpadding="0" cellspacing="6" width="100%" > 
      {1}
    </table
      {2}
      """.format('\n'.join(header), '\n'.join(table), '\n'.join(footer))
    return html  
    

  def row(self, rowdata, yellow='#FFD100', red='#FF1030'):
    '''
    creates a table row
    :type rowdata: list
    '''
    td = []
    bgcolor = ''
    try:
      if abs(float(rowdata[2])) > 2.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(yellow)
      if abs(float(rowdata[2])) > 3.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(red)
    except:
      pass
    for num, item in enumerate(rowdata): 
      try:
        # align right for numbers:
        float(item)
        td.append(r"""<td align='right' {0}> {1} </td>""".format(bgcolor, item))
      except:
        if item.startswith('-'):
          # only a minus sign
          td.append(r"""<td align='center'> {} </td>""".format(item))
        else:
          if num < 4:
            td.append(r"""<td align='right'> {} </td>""".format(item))
            continue
          # align left for words:
          td.append(r"""<td align='left'> {} </td>""".format(item))
    if not td:
      row = "<tr> No (disagreeable) restraints found in .lst file. </tr>"
    else:
      row = "<tr> {} </tr>".format(''.join(td))
    return row
    

ref = Refmod()  


if __name__ == '__main__':
  ref = Refmod()
  lst = ref.fileparser('d:\Programme\DSR\example\p21c-test.lst')
  htm = html_Table()
  header=['<h2>List of disagreeable restraints</h2>']
  footer = ['<br><br> Use "MORE 4" to get an extensive list of all restraints. ']
  print(htm.table_maker(header, lst, footer))