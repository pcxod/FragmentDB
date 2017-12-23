#FragmentDB Plugin

#Fragment Chooser
##Selection
Clicking on one of the list names selects a new fragment.

##Searching
Do an unsharp search with a partial name and **hit enter key**. 
For example a search for cf3 finds all the fragments containing cf3 in their 
name and similar.  

The search results in a shortened drop down menue.
  
##Reset
The reset button restores the full fragment list after a search.
  
##User Defined Database
The main database is write protected, but you can add and change as many 
fragments you like. They will be stored in a separate database in  
[Datadir()/db/user-fragment-database.sqlite](shell DataDir()/db). 
The user database will never be overwritten by plugin updates!  
Fragments created or changed by the user will get a \*user\* at the end of their names.  
    

#Picture
#Picture    
Displays a picture of the fragment.   
Click on the picture to display a magnified version of the picture.  
User defined fragments only have pictures if the user provided one.

#Properies
#Properies 
Define PART, free variables (FVAR) and the occupancy here. 
  
A free variable of 3 and an occupancy of 1 means 31.0 in SHELXL notation. 
-31.0 On the other hand is set with a free variable of -3. The corresponding free
variable is assigned automatically.  
All parameters will be displayed simultaneously through atomic labels.

#Properies 2
##Use a residue 
Disable or enable putting the fitted fragment into a residue. 
The residue class has to start with a letter and can be up to four  
characters. Usually, the class provided by the database is sufficient.

##Invert
Inverts the coordinates of the fragment. Useful to fit the inverted geometry  
of the fragment.

##DFIX
The supplied fragment normally has predefined restraints like 
SADI, FLAT and others.  
*DFIX* Generates DFIX/DANG restraints from the geometry of the fragment and 
replaces the predefined restraints. 


#Fit
##Fit!
The button <b>Fit!</b> starts the fit of the Fragment into the structure. 

##Edit
Edit the currenty selected or create new fragments.
  
###Change Picture
Select a picture for the fragment. Use at least 600dpi pictures.
  
###Drawstyle
Save a style file for chemdraw to draw molecular pictures that look 
the same as the predefined.
 
###Name
Name of the fragment

###Unit cell 
Cell parameters of the structure, where the fragment comes from. 
Use 1 1 1 90 90 90 for cartesian coordinates. The cell constants and the 
coordinates are converted into Cartesian coordinates during the adding to 
the database.

###Use selected atoms
Uses the selected atoms for the Atoms field and creates a picture from the 
fragment.

###Atoms 
insert atoms here like O2  8  1.3984  -0.3778   0.4922. You can also 
directly copy atoms from SHELX files or copy coordiinates from Avogadro. 
URL[http://avogadro.cc/, http://avogadro.cc]

###Restraints
Type any SHELX compatible restraint here (except SAME). 
The length of a line is unlimited.

###Residue
Define the residue class of the fragment (up to four characters, the first has to be a letter). 

##No Restraints
No restraints will be applied to the fitting fragment.

##Replace Mode
If enabled, all atoms inside PART 0 lying in a radius of 1.22 A around the 
atoms of the fitted fragment will be deleted. This is particular useful to 
replace atoms of a disordered structure just solved using SHELXT.  
It will not delete atoms in different PARTs than PART 0. 


#Results
#Results
This table shows the list of the most disagreeable restraints (in case there are some).  
A click on the respective line opens the editor window to adjust the restraints values.  
The long list shows ALL restraints and the short list only the most disargreeable ones.  
The error of a parameter should generally not be higher than 3sigma. A higher 
value means that the assumption of the target value originating from 
general observations is not met. The reason for this can be either, for 
example, a distorted molecular geometry like a bend aromatic ring or 
systematic errors in the data.  
Disagreeable restraints will be highlighted in yellow if the error is above 
2.5sigma. They will be highlighted in red if the error is above 3.5sigma. 
