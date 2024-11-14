#!/usr/bin/env python

########################################### Headers and Imports #################################################
import re
import os
import subprocess
import collections

from pymatgen.core import Composition, periodic_table
from pymatgen.entries.computed_entries import ComputedEntry

from pymatgen.analysis.phase_diagram import GrandPotPDEntry, PDEntry, GrandPotentialPhaseDiagram, PDPlotter, PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility


from mp_api.client import MPRester
from pymatgen.ext.matproj import MPRester

################################## Hydrogen and Oxygen Experimental Conditions ##################################
# Initialize compatibility module
compat = MaterialsProject2020Compatibility()

# Define dummy species for X
X = periodic_table.DummySpecie('X')

# Define compositions for gases
H2_Comp = Composition('H2')
O2_Comp = Composition('O2')
H2O_Comp = Composition('H2O')
CO_Comp = Composition('CO')
CO2_Comp = Composition('X')

# Define energies for hydrogen and oxygen under condition A
H_Ener_A = -4.024
O_Ener_A = -8.006 

# Define entries for gases under condition A
H2_Entry_A = ComputedEntry(H2_Comp, H_Ener_A * H2_Comp.num_atoms)
O2_Entry_A = ComputedEntry(O2_Comp, O_Ener_A * O2_Comp.num_atoms)
H2O_Entry_A = ComputedEntry(H2O_Comp, H_Ener_A * H2_Comp.num_atoms + O_Ener_A * 0.5 * O2_Comp.num_atoms)

# Define energies for hydrogen and oxygen under condition C
H_Ener_C = -4.997
O_Ener_C = -8.006

# Define entries for gases under condition C
H2_Entry_C = ComputedEntry(H2_Comp, H_Ener_C * H2_Comp.num_atoms)
O2_Entry_C = ComputedEntry(O2_Comp, O_Ener_C * O2_Comp.num_atoms)
H2O_Entry_C = ComputedEntry(H2O_Comp, H_Ener_C * H2_Comp.num_atoms + O_Ener_C * 0.5 * O2_Comp.num_atoms)

# Define entries for gases under condition X
O_Ener_X = -8.006
O2_Entry_X = ComputedEntry(O2_Comp, O_Ener_X * O2_Comp.num_atoms)
CO2_Ener_X = -25.556
CO_Ener_X = -20.232
CO_Entry_X = ComputedEntry(CO_Comp, CO_Ener_X)
CO2_Entry_X = ComputedEntry(CO2_Comp, CO2_Ener_X)

# Create lists of gas entries for different conditions
entriesGases_A = [H2_Entry_A, O2_Entry_A, H2O_Entry_A]
entriesGases_C = [H2_Entry_C, O2_Entry_C, H2O_Entry_C]
entriesGases_X = [O2_Entry_X, CO_Entry_X, CO2_Entry_X]

# Debugging: Print entries for debugging
# print("Entries for condition A:", entriesGases_A)
# print("----------------------------------------------------------------")
# print("Entries for condition C:", entriesGases_C)
# print("----------------------------------------------------------------")
# print("Entries for condition X:", entriesGases_X)

########################################### Material Entry ######################################################

# Define the compound you are working on
TestMat_Comp = Composition('Ba8Fe8O24')

# Debugging: get the atomic fraction of each element in the compound
# print(f"Reduced formula: {TestMat_Comp.reduced_formula}")
# print(f"Atomic fraction of Ba: {TestMat_Comp.get_atomic_fraction('Ba')}")
# print(f"Atomic fraction of Zr: {TestMat_Comp.get_atomic_fraction('Zr')}")
# print(f"Atomic fraction of O: {TestMat_Comp.get_atomic_fraction('O')}")
# print(f"Total number of atoms: {TestMat_Comp.num_atoms}")
# print(f"Elemental composition: {TestMat_Comp.get_el_amt_dict()}")

# Replace this with the computed energy for your material
TestMat_Ener = -268.28931146 # Placeholder, insert the correct value here

# Define computed entries for the material under condition A and C
TestMat_entry_A = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_A * 24)
TestMat_entry_C = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_C * 24)

# Debugging: print entries for debugging
# print("TestMat_entry_A:", TestMat_entry_A)
# print("TestMat_entry_C:", TestMat_entry_C)

# Create lists for VASP computed entries under different conditions
entries_VASP_A = [TestMat_entry_A]
entries_VASP_C = [TestMat_entry_C]

# Debugging: Print VASP computed entries
# print("entries_VASP_A:", entries_VASP_A)
# print("entries_VASP_C:", entries_VASP_C)

# Initialize MPRester to get material entries from the Materials Project
api = 'kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK'

# Automatically extract elements in TestMat_Comp
elements = [str(element) for element in TestMat_Comp.elements]

# Initialize MPRester to get material entries from the Materials Project
with MPRester(api) as api:
        entries_MP_Org_AC = api.get_entries_in_chemsys(
            elements + ['O', 'H']
        )
        entries_MP_Org_X = api.get_entries_in_chemsys(
            elements + ['O', 'C']
        )

# Debugging: Print fetched entries from Materials Project
# print("entries_MP_Org_AC:", entries_MP_Org_AC)
# print("entries_MP_Org_X:", entries_MP_Org_X)

# Process the entries using the compatibility module
entries_MP_Org_AC = compat.process_entries(entries_MP_Org_AC)
entries_MP_Org_X = compat.process_entries(entries_MP_Org_X)

# Debugging: Print processed entries
# print("Processed entries_MP_Org_AC:", entries_MP_Org_AC)
# print("Processed entries_MP_Org_X:", entries_MP_Org_X)

# Combine the entries with VASP computed entries
all_entries_A = entries_MP_Org_AC + entries_VASP_A
all_entries_C = entries_MP_Org_AC + entries_VASP_C
entriesTotal_X = entries_MP_Org_X

# Debugging: Print reduced formulas for condition A entries horizontally
# print("Entries for condition A:")
# print(" | ".join([entry.composition.reduced_formula for entry in all_entries_A]))
##################################### Eliminate H2 and O2 MP Entries ###########################################

# Filter out H2, O2, and H2O from the list of all entries
H2_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'H2']
O2_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'O2']
H2O_entries_A = [e for e in all_entries_A if e.composition.reduced_formula == 'H2O']

# Debugging: Print filtered entries for condition A
# print("Filtered H2 entries for condition A:", H2_entries_A)
# print("Filtered O2 entries for condition A:", O2_entries_A)
# print("Filtered H2O entries for condition A:", H2O_entries_A)

# Filter out H2, O2, and H2O from the list of all entries
H2_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'H2']
O2_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'O2']
H2O_entries_C = [e for e in all_entries_C if e.composition.reduced_formula == 'H2O']

# Debugging: Print filtered entries for condition C
# print("Filtered H2 entries for condition C:", H2_entries_C)
# print("Filtered O2 entries for condition C:", O2_entries_C)
# print("Filtered H2O entries for condition C:", H2O_entries_C)

# Filter out H2, O2, and H2O from the list of all entries
eliminate_AC = ['H2', 'O2', 'H2O']

all_entries_A = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_A))
all_entries_A =all_entries_A +entriesGases_A

all_entries_C = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC,  all_entries_C))
all_entries_C =all_entries_C +entriesGases_C

# Debugging: check if the entries contain the gases 
# for entry in entriesGases_A:
#     print(entry.composition.reduced_formula)
# print("************************************************************************************************")

# Debugging: Print entries to ensure they contain oxygen
# print("Entries for condition A:")
# print(" | ".join([entry.composition.reduced_formula for entry in all_entries_A]))

# Debugging: Confirm the presence of BaZrO3 in the list of entries
# found_BaZrO3 = False
# for entry in all_entries_A:
#     if entry.composition.reduced_formula == 'BaZrO3':
#         print("Found BaZrO3:", entry)
#         found_BaZrO3 = True
# if not found_BaZrO3:
#     print("BaZrO3 not found in entries for condition A")

# Debugging: Print entries to ensure they contain oxygen
# print("Entries for condition C:")
# print(" | ".join([entry.composition.reduced_formula for entry in all_entries_C]))

########################## Calculate Convex Hull under Environmental Conditions A+C #############################

# Locked chemical potentials for condition A
locked_Chem_Potential_A = {'H2': H_Ener_A*2,'O2': O_Ener_A*2}
# print("Chemical potential for condition A:", locked_Chem_Potential_A)
# print("entriesGases_A: ", " | ".join(str(entry) for entry in entriesGases_A))
# Debugging: Output locked chemical potential for condition A
# print("Chemical potential for condition A:", locked_Chem_Potential_A)

# Debugging: Print entries to ensure they contain oxygen
# for entry in all_entries_A:
#     if 'O' in entry.composition.reduced_formula:
#         print("O entry found:", entry.composition.reduced_formula)

# Debugging: Print entries to ensure they contain oxygen
# print("Entries for condition A before creating phase diagram:")
# for entry in all_entries_A:
#     print(entry.composition.reduced_formula)

# Create phase diagram for condition A by using GrandPotentialPhaseDiagram
pd_A = GrandPotentialPhaseDiagram(all_entries_A, locked_Chem_Potential_A)

# Debugging: Print entries to ensure they contain oxygen
# print("Entries for condition A after creating phase diagram:", pd_A.all_entries)

# Output phase diagram and convex hull energy for condition A
print('************************************************************************************************')
print("Phase diagram for Hydrogen-rich condition:", pd_A)
print("every energy is per atom: ",TestMat_entry_A.energy / 16)
print(f"Hull energy for compound Ba8Zr8O24: {pd_A.get_hull_energy(TestMat_entry_A.composition)/16}")
print("Energy above convex hull:",TestMat_entry_A.energy / 16 - pd_A.get_hull_energy_per_atom(TestMat_entry_A.composition))
# plotter = PDPlotter(pd_A)
# plotter.show()

# Write results to file
with open('Anode_Ba8Zr8O24.txt', 'w') as out_file:
    print(TestMat_Comp, file=out_file)
    print(TestMat_entry_A.energy / 16 - pd_A.get_hull_energy_per_atom(TestMat_entry_A.composition), file=out_file)
    print('************************************************************************************************')
# Locked chemical potentials for condition C
locked_Chem_Potential_C = {'O2': O_Ener_C*2,'H2': H_Ener_C*2}
# print("Chemical potential for condition C:", locked_Chem_Potential_C)
# print("entriesGases_C: ", " | ".join(str(entry) for entry in entriesGases_C))
pd_C = GrandPotentialPhaseDiagram(all_entries_C, locked_Chem_Potential_C)
    
# Output phase diagram and convex hull energy for condition C
print("Phase diagram for Oxygen-rich condition:",pd_C)
print("every energy is per atom: ", TestMat_entry_C.energy / 16)
print(f"Hull energy for compound Ba8Zr8O24: {pd_C.get_hull_energy(TestMat_entry_C.composition)/16}")
print("Energy above convex hull:",TestMat_entry_C.energy / 16 - pd_C.get_hull_energy_per_atom(TestMat_entry_C.composition))
# plotter = PDPlotter(pd_C)
# plotter.show()
# Write results to file
with open('Cathode_Ba8Zr8O24.txt', 'w') as out_file:
    print(TestMat_Comp, file=out_file)
    print(TestMat_entry_C.energy / 16 - pd_C.get_hull_energy_per_atom(TestMat_entry_C.composition), file=out_file)
    print('************************************************************************************************')
########################################## Filtering CO2 as Element X ###########################################

# Modify composition to account for CO2 (as X)
for j in range(0, len(entriesTotal_X)):
    CurrentEntry = entriesTotal_X[j]
    get_Composition = CurrentEntry.composition

    # Modify composition text
    with open('composition.txt', 'w') as out_file:
        originalLine = str(get_Composition)
        fixedLine1 = originalLine.replace('\n', 'Ba8Zr8O24')
        fixedLine2 = re.sub("[A-Za-z]+", lambda ele: " " + ele[0] + " ", fixedLine1)
        print(fixedLine2, file=out_file)

    # Execute external script to split composition
    subprocess.run('./split_composition.sh')

    with open('./num_Lines.txt', 'r') as out_file:
        numLinesO = out_file.readlines()
        numLines = numLinesO[0]
        # numLines = numLines.replace('\n', 'Ba8Zr8O24')
        numLines = int(numLines)

    numO = 0
    numC = 0
    holder = []
    holder = [0 for i in range(numLines)]

    with open('composition.txt','r') as out_file:
        holderOut = out_file.readlines()
        for k in range(0, numLines):
            holder[k] = holderOut[k]
            holder[k] = holder[k].replace('\n', 'Ba8Zr8O24')

            if holder[k] == 'C':
                numCs = holderOut[k + 1]
                numCs = str(numCs.replace('\n', 'Ba8Zr8O24'))
                numCs = int(numCs)
            elif holder[k] == 'O':
                numOs = holderOut[k + 1]
                numOs = str(numOs.replace('\n', 'Ba8Zr8O24'))
                numOs = int(numOs)

    C_exist_check = 'numC' in locals()
    O_exist_check = 'numO' in locals()

    if C_exist_check == False:
        numC = int(0)

    if O_exist_check == False:
        numO = int(0)

    if numC == 0 and numO == 0:
        entriesTotal_X[j] = CurrentEntry
    elif numC > 0 and numO == 0:
        entriesTotal_X[j] = CurrentEntry
    elif numC == 0 and numO > 0:
        entriesTotal_X[j] = CurrentEntry
    elif numC > 0 and numO < 2:
        entriesTotal_X[j] = CurrentEntry
    elif numC > 0 and numO > 0:
        numX = 0
        Mod_Comp = 'Ba'
        if(2*numC == numO):
            numX = numC
            g=1
            while g< numLines+1:
                if holder[g-1]=='C':
                    Mod_Comp = Mod_Comp + 'X' + str(numX)
                    g=g+2
                elif holder[g-1]=='O':
                    g=g+2
                else:
                    Mod_Comp = holder[g-1]+holder[g]+Mod_Comp
                    g=g+2
        elif(2*numC > numO):
            while True:
                numO = numO - 2
                if numO < 0:
                    numO = numO + 2
                    break
                numX = numX + 1
                numC = numC - 1
            g=1
            while g< numLines+1:
                if holder[g-1]=='C':
                    Mod_Comp = Mod_Comp + 'X' + str(numX)
                    Mod_Comp = Mod_Comp +'C'+str(numC)
                    g=g+2
                elif holder[g-1]=='O':
                    if numO == 0:
                        g=g+2
                    elif numO > 0:
                        Mod_Comp = Mod_Comp + 'O' + str(numO)
                        g=g+2
                else:
                    Mod_Comps=holder[g-1]+holder[g]+Mod_Comp
                    g=g+2
        elif(2*numC < numO):    
            numX = numC
            numO = numO - 2*numC
            g=1
            while g<numLines+1:
                if holder[g-1]=='C':
                    Mod_Comp = Mod_Comp + 'X' 
                    Mod_Comp = Mod_Comp + str(numC)
                    g=g+2
                elif holder[g-1]=='O':
                    Mod_Comp = Mod_Comp + 'O'
                    Mod_Comp = Mod_Comp + str(numO)
                    g=g+2
                else:
                    Mod_Comp = holder[g-1]+holder[g]+Mod_Comp
                    g=g+2
            energyEntry = CurrentEntry.energy
            make_Mod_Entry = ComputedEntry(Mod_Comp, energyEntry)
            entriesTotal_X[j] = make_Mod_Entry

################################### CO(Element Z) Checks ###################################################
for j in range(0,len(entriesTotal_X)):

    CurrentEntry = entriesTotal_X[j]
    get_Composition = CurrentEntry.composition

    with open('composition.txt', 'w') as out_file:
        originalLine = str(get_Composition)
        fixedLine1 = originalLine.replace('\n', 'Ba8Zr8O24')
        fixedLine2 = re.sub("[A-Za-z]+", lambda ele: " " + ele[0] + " ", fixedLine1)
        print(fixedLine2, file=out_file)

    subprocess.run('./split_composition.sh')

    with open('./num_Lines.txt', 'r') as out_file:
        numLinesO = out_file.readlines()
        numLines = numLinesO[0]
        # numLines = numLines.replace('\n', 'Ba8Zr8O24')
        numLines = int(numLines)

    numO = 0
    numC = 0
    holder = []
    holder = [0 for chi in range(numLines)]

    with open('composition.txt','r') as out_file:
        holderOut = out_file.readlines()

        for k in range(0, numLines):
            holder[k] = holderOut[k]
            holder[k] = holder[k].replace('\n', 'Ba8Zr8O24')

            if holder[k] == 'C':
                numCs = holderOut[k + 1]
                numCs = str(numCs.replace('\n', 'Ba8Zr8O24'))
                numCs = int(numCs)
            elif holder[k] == 'O':
                numOs = holderOut[k + 1]
                numOs = str(numOs.replace('\n', 'Ba8Zr8O24'))
                numOs = int(numOs)
        
        C_exist_check = 'numC' in locals()
        O_exist_check = 'numO' in locals()

        if C_exist_check == False:
            numC = int(0)

        if O_exist_check == False:
            numO = int(0)
        
        if numC>0 and numO==0:
            numZ = 0
            Mod_Comp = 'Ba8Zr8O24'
            if(numC == numO):
                numZ = numC
                g=1
                while g< numLines+1:
                    if holder[g-1]=='C':
                        Mod_Comp = Mod_Comp + 'Z' + str(numZ)
                        g=g+2
                    elif holder[g-1]=='O':
                        g=g+2
                    else:
                        Mod_Comp = holder[g-1]+holder[g]+Mod_Comp
                        g=g+2
            elif(numC > numO):
                while True:
                    numC = numC - 1
                    if numC < 0:
                        numC = numC + 1
                        break
                    numZ = numZ + 1
                    numO = numO - 1
                g=1
                while g< numLines+1:
                    if holder[g-1]=='C':
                        Mod_Comp = Mod_Comp + 'Z' + str(numZ)
                        Mod_Comp = Mod_Comp +'C'+str(numC)
                        g=g+2
                    elif holder[g-1]=='O':
                        if numO == 0:
                            g=g+2
                        elif numO > 0:
                            Mod_Comp = Mod_Comp + 'O' + str(numO)
                            g=g+2
                    else:
                        Mod_Comps=holder[g-1]+holder[g]+Mod_Comp
                        g=g+2
                
                energyEntry = CurrentEntry.energy
                make_Mod_Entry = ComputedEntry(Mod_Comp, energyEntry)
                entriesTotal_X[j] = make_Mod_Entry
                print(entriesTotal_X[j])

#################################### Gas Phase Correction CO2 ###################################################

CO_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'CO']
CO2_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'X']
O2_entries_X = [e for e in entriesTotal_X if e.composition.reduced_formula == 'O2']

eliminate_X = ['CO', 'X', 'O2']
all_entries_X = list(filter(lambda e: e.composition.reduced_formula not in eliminate_X, entriesTotal_X))

TestMat_entry_X = ComputedEntry(TestMat_Comp, TestMat_Ener - O_Ener_X * 24)
entries_VASP_X = [TestMat_entry_X]
all_entries_X = all_entries_X + entries_VASP_X + entriesGases_X

############################ Calculate Convex Hull under Environmental Conditions ###############################

# Locked chemical potentials for condition X
locked_Chem_Potential_X = {'X': CO2_Ener_X,'O2': O_Ener_X*2}
# print("Chemical potential for condition X:", locked_Chem_Potential_X)
# print("entriesGases_X: ", " | ".join(str(entry) for entry in entriesGases_X))
pd_X = GrandPotentialPhaseDiagram(all_entries_X, locked_Chem_Potential_X)

# Output phase diagram and convex hull energy for condition X
print("Phase diagram for CO2-rich condition:",pd_X)
print("every energy is per atom: ",TestMat_entry_X.energy / 16)
print(f"Hull energy for compound Ba8Zr8O24: {pd_X.get_hull_energy(TestMat_entry_X.composition)/16}")
print("Energy above convex hull:",TestMat_entry_X.energy / 16 - pd_X.get_hull_energy_per_atom(TestMat_entry_X.composition))
# plotter = PDPlotter(pd_X)
# plotter.show()
# Write results to file
with open('CO2_Ba8Zr8O24_X.txt', 'w') as out_file:
    print(TestMat_Comp, file=out_file)
    print(TestMat_entry_X.energy / 16 - pd_X.get_hull_energy_per_atom(TestMat_entry_X.composition), file=out_file)
print('************************************************************************************************')