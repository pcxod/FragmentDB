
'''
handles the fragment atoms of a new/edited fragment
C1 1 x y z 11.0 U U =   U U U U
C1 1 x y z 11.0
next atom after the fifth list element has to start with a character
and not "=", otherwise it is a U value


      for num, i in enumerate(atoms):
    if num >= 5 and not i[0].isalpha():
        del atoms[6]
        continue
      atomline = atoms[i:i+5]
      if len(atomline) < 5:
        #print("Parameter of Atom {} missing".format(atomline[0]))
        continue
      atlines.append(' '.join(atomline))

take the firat 5 elements from the string, delete them, check ich first char from
the next list element isalpha(), otherwise delete the next element as long as the 
first char of the next string isalpha().
'''
    
atoms = '''C4    1    0.282212   0.368636    0.575493   
C7    4    0.358324    0.417816    0.593763    
F8    4    0.351897    0.315202    0.582798    11.00000    0.03627    0.03415 =
         0.04126    0.00790   -0.02038    0.00477
GA1   6    0.639512    0.561734    0.237759    11.00000    0.02407    0.02504 =
         0.02490    0.00000   -0.00156    0.00159
AL1   5    0.064277    0.260189    0.234
'''

atlines = []
atoms=atoms.split()
atline = []
for num, i in enumerate(atoms):
  atline.append(i)
  try:
    # cut the long list in chuncs of atoms:
    if num > 1 and atoms[num+1][0].isalpha():
      atlines.append(atline)
      atline = []
  except IndexError:
    # the last atom has no num+1
    atlines.append(atline)
# go through all atoms and cut their line to 5 list elements At SFAC x y z:
for num, line in enumerate(atlines):
  if len(line) > 5:
    atlines[num] = line[:5]
  if len(line) < 5:
    # too short, parameters missing
    print('Bad atom line found!! Parameter(s) missing.')
  for x in line[1:5]:
    # check if each is a rea number except for the atom:
    try:
      float(x)
    except:
      # if something is wrong, determine the bad guy:
      for i in x:
        if not i.isdigit() and i != '.':
          print('Bad charachter {} in line.'.format(i))
          continue
      print('Bad atom line found!')
# return atlines

for i in atlines:
  print(i)
  
  