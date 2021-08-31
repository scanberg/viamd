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

## Animation ##
- Show time and frames (total number)

## Menu and Windows ##
- Add Help
- Add Hierarchical info of dataset (chains, residues, atoms, element types what ever)

## Navigation ##
- Animate rotation
- Smooth rotation overall

## Properties ##
- rdf

- Fix GUI, the god damn columns are messed up
- Fix bug with range based selection
- <span style="color:red">Visualization (Show cloud)</span>
- <span style="color:red">Coloring (Color cloud based on property value)</span>

## Timeline ##
- Add units on hover
- Export to ascii table

EXAMPLE OF EXPORT:
[//]: #  asd: distance resatom(resname(ALA), 6) resatom(resname(ALA), 4)
[//]: # asd1: distance resatom(residue(1), 6) resatom(residue(1), 4)
[//]: # asd2: distance resatom(residue(2), 6) resatom(residue(2), 4)
[//]: #   a1: angle atom(73) atom(75) atom(81)

time    asd     asd.1     asd.2       a1


## Distributions ##
- CTRL + CLICK AND DRAG = SELECT RANGE
- Periodic Selection
- Value??? Assert this for all types
- Add Units
- Export to table data

## Volume ##
Intersection of filtered properties

## Hydrogen Bonds ##
Fix computation assert validity
Compute per shown frame

## Ramachandran ##
- Statistics (How many in each bin)
- Make aggregate distribution for selection / highlight?
- Follow color scheme (Blue is selected) (Yellow is highlight)


## Selection ##
- Fix selection of hidden particles (However, do show hidden selections??)