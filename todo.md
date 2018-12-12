# TODO VIAMD #

## Loading Structures ##

### PDB ###
- CG
    - Let user specify radii for beads
- List of different residues
- Supply time between frames

### GRO ###
- CG
    - Let user specify radii for beads
- List different residues
- Load Trajectory Button
- Supply time between frames

## Trajectory Explorer ##
- Change label of window
- Show time and frames (total number)
- If stopped at last frame, click on play to start from beginning

## Menu and Windows ##
- Add Help
- Add Hierarchical info of dataset (chains, residues, atoms, element types what ever)

## Navigation ##
- Double click -> Center on atom
- Use animation

## Representation ##
- Fix Bugs
- Range should be :, Examples: 1:10, 10:*, *
- Filter
    - Backbone [ ]
    - Water    [x]

## Properties ##
- distance
- angle
- dihedral
- rdf
- rmsd

- com
- resatom(residues, idx)
- atom(int)
- resname(string)

- Fix bug with range based selection
- <span style="color:red">Visualization (Show cloud)</span>
- <span style="color:red">Coloring (Color cloud based on property value)</span>

## Timeline ##
- Add units on hover
- Export to ascii table

[//]: #  asd: distance resatom(resname(ALA), 6) resatom(resname(ALA), 4)
[//]: # asd1: distance resatom(residue(1), 6) resatom(residue(1), 4)
[//]: # asd2: distance resatom(residue(2), 6) resatom(residue(2), 4)
[//]: #   a1: angle atom(73) atom(75) atom(81)

time    asd     asd.1     asd.2       a1


## Distributions ##
- Periodicity Selection
- Value??? Assert this for all types
- Add Units
- Export to ascii table


## Volume ##
Intersection of filtered properties

## Hydrogen Bonds ##
Fix computation assert validity
Compute per shown frame

## Ramachandran ##
- Fix size of image
- Show values of X and Y axis
- Show phi and psi on axis
- Filter (Same as representation)
- Make nicer over all
- Statistics (How many in each bin)