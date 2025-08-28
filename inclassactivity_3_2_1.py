##########################################################
## BME 494/598:  In-class activity 3.2 #1               ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed for: BME 494/598 Applied Programming      ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

## Your name
# <Rori hoover>


## Import required packages
import json


## 1. Load up periodic table of elements data
# Deliverables:
#   a. Open the 'periodic_table_lookup.json' file
#   b. Load periodic table of elements data into variable called 'ptoe'
with open('/Users/rorihoover/Downloads/periodic_table_lookup.json', 'r') as file:
    ptoe_all = json.load(file)

#   c. Print the keys for ptoe_all
print(ptoe_all.keys())

#   d. Print the key 'order'
print(ptoe_all["order"])

#   e. Print the key 'carbon'
print("QUESTION E")
print(ptoe_all["carbon"])

#   f. Print the atomic symbol for 'carbon'
print("Atomic symbol of carbon")
print(ptoe_all["carbon"]["symbol"])

#   g. Print the atomic mass for 'carbon'
print("atomic mass of carbon", ptoe_all["carbon"]["atomic_mass"] )


## 2. (5pts) Build a dictionary that is keyed off atomic symobl and has atomic mass as a value
#     Molecules will be given as dictionaries with atomic symbol as key and value being the
#     number of atoms. The ptoe_all is keyed off the full atom name, and we want to have the
#     atomic mass keyed off the atomic symbol. {'C': 12.011, ... }
# Deliverables:
#   a. Initialize a new dictionary 'ptoe'
#   b. Iterate through all 'ptoe_all' atoms (easiest to use ptoe_all['order'] as the iterable)
#   c. Add each atom into 'ptoe' with atomic symbol ('symbol') as the key and atomic mass
#      ('atomic_mass') as the value

ptoe = {}
for i in ptoe_all['order']:
    symbol = ptoe_all[i]["symbol"]
    mass = ptoe_all[i]["atomic_mass"]
    ptoe[symbol] = mass

print(ptoe)
#   d. Print the atomic mass of 'C' using 'ptoe', and check if equals 12.011
#   e. Print the atomic mass of 'Fe' using 'ptoe', and check if equals 55.8452
#   f. Print the atomic mass of 'U' using 'ptoe', and check if equals 238.028913
print(ptoe["C"])
print(ptoe["C"]== 12.011)
print(ptoe["Fe"])
print(ptoe["Fe"]== 55.8452)
print(ptoe["U"])
print(ptoe["U"]==238.028913)

## 3. (5pts) Define a function 'molecular_weight' to find the molecular weight for molecules
# Deliverables:
#   a. Write a function 'molecular_weight' that:
#      - Takes two arguments 'molecule' and 'ptoe'
#        a. 'molecule' argument is a dictionary of atoms in the molecule with the value being the number of atoms
#        b. 'ptoe' the dictionary you developed above with atomic symbol as keys and atomic mass as values
#      - Use the atomic weights and the molecular composition to compute the molecular weight
#      - Return the molecular weight
#      - Docstrings to document the function
def molecular_weight(molecule, ptoe):
    """
        Function calculates the molecular weight of molecules 
    Parameters:
        moleclue: a dictionary of atoms in the molecule with the value being the number of atoms
        ptoe: dictionary with atomic symbol as keys and atomic mass as values
    Returns:
        the molecular weight
    """
    result = 0.0
    for i in molecule.keys():
        molecule_count = molecule[i]
        mass = ptoe[i] * molecule_count
        result = result + mass
    return result

Glycine = {"C": 2, "H":5,"N": 1, "O": 2}
Glucose = {"C":6,"H":12,"O":6}
Palmitic_acid = {"C":16,"H":32,"O":2}
ATP = {"C": 10, "H": 16, "N": 5, "O": 13, "P": 3}
Dichlorodifluoromethane = {"C": 1, "Cl": 2, "F": 2}
Selenocysteine = {"C": 3, "H": 7, "N": 1, "O": 2, "Se": 1}
Heme_B = {"C": 34, "H": 32, "O": 4, "N": 4, "Fe": 1}


print("Glycine:", molecular_weight(molecule = Glycine, ptoe = ptoe ))
print("Glucose:", molecular_weight(molecule = Glucose, ptoe = ptoe ))
print("Palmitic acid:", molecular_weight(molecule = Palmitic_acid, ptoe = ptoe ))
print("ATP:", molecular_weight(molecule = ATP, ptoe = ptoe ))
print("Dichlorodifluoromethane:", molecular_weight(molecule = Dichlorodifluoromethane, ptoe = ptoe ))
print("Selenocysteine:", molecular_weight(molecule = Selenocysteine, ptoe = ptoe ))
print("Heme B:", molecular_weight(molecule = Heme_B, ptoe = ptoe ))






