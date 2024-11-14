#!/usr/bin/env python

# Headers and Imports
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from mp_api.client import MPRester
import sys
import os


# Global Initialization Hydrogen and Oxygen Experimental Conditions
def initialize_global_variables():
    """
    Initializes global gas entries and chemical potentials for conditions A, C, and X.
    This function sets up the required parameters for each environment (A, C, and X) by
    defining energies, compositions, and entries for relevant gases.
    """

    global H_Ener_A, O_Ener_A, H_Ener_C, O_Ener_C, O_Ener_X, CO2_Ener_X, CO_Ener_X
    global H2_Comp, O2_Comp, H2O_Comp, CO_Comp, CO2_Comp
    global entriesGases_A, entriesGases_C, entriesGases_X
    global locked_Chem_Potential_A, locked_Chem_Potential_C, locked_Chem_Potential_X

    # Define energies for different environments
    H_Ener_A = -4.024
    O_Ener_A = -8.006
    H_Ener_C = -4.997
    O_Ener_C = -6.166
    O_Ener_X = -6.320
    CO2_Ener_X = -25.556
    CO_Ener_X = -20.232

    # Define compositions for gases
    H2_Comp = Composition("H2")
    O2_Comp = Composition("O2")
    H2O_Comp = Composition("H2O")
    CO_Comp = Composition("CO")
    CO2_Comp = Composition("X")

    # Define entries for each condition
    entriesGases_A = [
        ComputedEntry(H2_Comp, H_Ener_A * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_Ener_A * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_Ener_A * H2O_Comp.num_atoms + O_Ener_A * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_C = [
        ComputedEntry(H2_Comp, H_Ener_C * H2_Comp.num_atoms),
        ComputedEntry(O2_Comp, O_Ener_C * O2_Comp.num_atoms),
        ComputedEntry(H2O_Comp, H_Ener_C * H2O_Comp.num_atoms + O_Ener_C * 0.5 * O2_Comp.num_atoms)
    ]
    entriesGases_X = [
        ComputedEntry(O2_Comp, O_Ener_X * O2_Comp.num_atoms),
        ComputedEntry(CO_Comp, CO_Ener_X),
        ComputedEntry(CO2_Comp, CO2_Ener_X)
    ]

    # Define chemical potentials
    locked_Chem_Potential_A = {'H2': H_Ener_A * 2, 'O2': O_Ener_A * 2}
    locked_Chem_Potential_C = {'H2': H_Ener_C * 2, 'O2': O_Ener_C * 2}
    locked_Chem_Potential_X = {'O2': O_Ener_X * 2, 'X': CO2_Ener_X}
    # locked_Chem_Potential_A = {'H': H_Ener_A , 'O': O_Ener_A }
    # locked_Chem_Potential_C = {'H': H_Ener_C , 'O': O_Ener_C }
    # locked_Chem_Potential_X = {'O': O_Ener_X , 'X': CO2_Ener_X}


# Material Entry
def prepare_material_entries(api, TestMat_Comp, TestMat_Ener):
    """
    Prepare material entries and fetch entries from the Materials Project for a given material.

    Args:
        api: API key for the Materials Project.
        TestMat_Comp (str): Composition of the material.
        TestMat_Ener (float): Energy of the material.

    Returns:
        tuple: List of entries for Condition A, list of entries for Condition C, list of entries for Condition X,
               ComputedEntry object for the material under Condition A,
               ComputedEntry object for the material under Condition C,
               ComputedEntry object for the material under Condition X.
    """

    # Define the material's composition
    TestMat_Comp = Composition(TestMat_Comp)

    # Define computed entries for the material under conditions A, C and X with unique entry_ids
    TestMat_entry_A = ComputedEntry(TestMat_Comp, TestMat_Ener)
    TestMat_entry_C = ComputedEntry(TestMat_Comp, TestMat_Ener)
    TestMat_entry_X = ComputedEntry(TestMat_Comp, TestMat_Ener)

    # Initialize compatibility module
    compat = MaterialsProject2020Compatibility()

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

    # Process entries with the compatibility module
    entries_MP_Org_AC = compat.process_entries(entries_MP_Org_AC)
    entries_MP_Org_X = compat.process_entries(entries_MP_Org_X)

    # Combine the fetched entries with the VASP computed entries for Conditions A and C
    entries_VASP_A = [TestMat_entry_A]
    entries_VASP_C = [TestMat_entry_C]

    all_entries_A = entries_MP_Org_AC + entries_VASP_A
    all_entries_C = entries_MP_Org_AC + entries_VASP_C
    entriesTotal_X = entries_MP_Org_X

    return all_entries_A, all_entries_C, entriesTotal_X, TestMat_entry_A, TestMat_entry_C, TestMat_entry_X


# Condition A
def calculate_phase_diagram_condition_A(all_entries_A, entriesGases_A, locked_Chem_Potential_A, TestMat_entry_A):
    """
    Calculate the phase diagram and energy above the hull for Condition A (Hydrogen-rich environment).

    Args:
        all_entries_A (list): List of ComputedEntry objects for Condition A.
        entriesGases_A (list): List of ComputedEntry objects for gases under Condition A.
        locked_Chem_Potential_A (dict): Locked chemical potentials for Condition A.
        TestMat_entry_A (ComputedEntry): ComputedEntry object for the material under Condition A.

    Returns:
        tuple: Phase diagram for Condition A, energy per atom, formation energy, and energy above hull.
    """

    # Define species to eliminate
    eliminate_AC = ['H2', 'O2', 'H2O']

    # Filter out these species from all_entries_A
    all_entries_A = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_A))
    all_entries_A = all_entries_A + entriesGases_A
    
    # Create phase diagram for Condition A
    pd_A = GrandPotentialPhaseDiagram(all_entries_A, locked_Chem_Potential_A)

    # Calculate energy per atom
    energy_per_atom_A = TestMat_entry_A.energy / 16

    # Get the grand potential entry
    gpe = next((e for e in pd_A.all_entries if e.original_entry == TestMat_entry_A), None)

    # Calculate formation energy
    formation_energy_A = pd_A.get_form_energy_per_atom(gpe)

    # # Calculate energy above hull
    energy_above_hull_A = pd_A.get_e_above_hull(gpe)
    # energy_above_hull_A = TestMat_entry_A.energy / 16 - pd_A.get_hull_energy_per_atom(TestMat_entry_A.composition)

    return pd_A, energy_per_atom_A, formation_energy_A, energy_above_hull_A


# Condition C
def calculate_phase_diagram_condition_C(all_entries_C, entriesGases_C, locked_Chem_Potential_C, TestMat_entry_C):
    """
    Calculate the phase diagram and energy above the hull for Condition C (Oxygen-rich environment).

    Args:
        all_entries_C (list): List of ComputedEntry objects for Condition C.
        entriesGases_C (list): List of ComputedEntry objects for gases under Condition C.
        locked_Chem_Potential_C (dict): Locked chemical potentials for Condition C.
        TestMat_entry_C (ComputedEntry): ComputedEntry object for the material under Condition C.

    Returns:
        tuple: Phase diagram for Condition C, energy per atom, formation energy, and energy above hull.
    """

    # Define species to eliminate
    eliminate_AC = ['H2', 'O2', 'H2O']

    # Filter out these species from all_entries_C
    all_entries_C = list(filter(lambda e: e.composition.reduced_formula not in eliminate_AC, all_entries_C))
    all_entries_C = all_entries_C + entriesGases_C

    # Create phase diagram for Condition C
    pd_C = GrandPotentialPhaseDiagram(all_entries_C, locked_Chem_Potential_C)

    # Calculate energy per atom
    energy_per_atom_C = TestMat_entry_C.energy / 16

    # Get the grand potential entry
    gpe = next((e for e in pd_C.all_entries if e.original_entry == TestMat_entry_C), None)

    # Calculate formation energy
    formation_energy_C = pd_C.get_form_energy_per_atom(gpe)

    # Calculate energy above hull
    energy_above_hull_C = pd_C.get_e_above_hull(gpe)
    #energy_above_hull_C = TestMat_entry_C.energy / 16 - pd_C.get_hull_energy_per_atom(TestMat_entry_C.composition)
    
    return pd_C, energy_per_atom_C, formation_energy_C, energy_above_hull_C


# Condition X
def calculate_phase_diagram_condition_X(entriesTotal_X, entriesGases_X, locked_Chem_Potential_X, TestMat_entry_X):
    """
    Calculate the phase diagram and energy above the hull for Condition X (CO2-rich environment).

    Args:
        entriesTotal_X (list): List of ComputedEntry objects for Condition X.
        entriesGases_X (list): List of ComputedEntry objects for gases under Condition X.
        locked_Chem_Potential_X (dict): Locked chemical potentials for Condition X.
        TestMat_entry_X (ComputedEntry): ComputedEntry object for the material under Condition X.

    Returns:
        tuple: Phase diagram for Condition X, energy per atom, formation energy, and energy above hull.
    """

    # Define species to eliminate
    eliminate_X = ['CO', 'X', 'O2']

    # Filter out these species from entriesTotal_X
    all_entries_X = list(filter(lambda e: e.composition.reduced_formula not in eliminate_X, entriesTotal_X))

    # Add back entries for Condition X
    entries_VASP_X = [TestMat_entry_X]
    all_entries_X = all_entries_X + entries_VASP_X + entriesGases_X

    # Create phase diagram for Condition X
    pd_X = GrandPotentialPhaseDiagram(all_entries_X, locked_Chem_Potential_X)

    # Calculate energy per atom
    energy_per_atom_X = TestMat_entry_X.energy / 16

    # Get the grand potential entry
    gpe = next((e for e in pd_X.all_entries if e.original_entry == TestMat_entry_X), None)

    # Calculate formation energy
    formation_energy_X = pd_X.get_form_energy_per_atom(gpe)

    # Calculate energy above hull
    energy_above_hull_X = pd_X.get_e_above_hull(gpe)
    # energy_above_hull_X = TestMat_entry_X.energy / 16 - pd_X.get_hull_energy_per_atom(TestMat_entry_X.composition)

    return pd_X, energy_per_atom_X, formation_energy_X, energy_above_hull_X


# Main Function
def main():
    # Determine input compositions
    if len(sys.argv) < 2:
        # No command-line arguments provided; attempt to read from 'composition.txt'
        input_file = 'composition.txt'
        with open(input_file, 'r') as f:
            compositions = [line.strip() for line in f if line.strip()]
    else:
        # Compositions provided as command-line arguments
        compositions = sys.argv[1:]

    # Initialize gas entries and chemical potentials
    initialize_global_variables()

    # API key for Materials Project
    api_key = 'kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK' 

    # Iterate over each composition
    for comp_str in compositions:
        # Define the composition
        TestMat_Comp = comp_str

        # Placeholder for TestMat_Ener; replace with actual energy value or a method to obtain it
        TestMat_Ener = -334.17059441 # Ba8Zr8O24
        # TestMat_Ener = -251.41021058 # Ba8Co8O24
        # TestMat_Ener = -268.01045909 # Ba8Fe8O24

        print(f"\nProcessing composition: {comp_str}")

        # Prepare material entries
        all_entries_A, all_entries_C, entriesTotal_X, TestMat_entry_A, TestMat_entry_C, TestMat_entry_X = prepare_material_entries(
            api=api_key,
            TestMat_Comp=TestMat_Comp,
            TestMat_Ener=TestMat_Ener
        )

        # Calculate phase diagram for Condition A
        pd_A, energy_per_atom_A, formation_energy_A, energy_above_hull_A = calculate_phase_diagram_condition_A(
            all_entries_A,
            entriesGases_A,
            locked_Chem_Potential_A,
            TestMat_entry_A
        )

        # Calculate phase diagram for Condition C
        pd_C, energy_per_atom_C, formation_energy_C, energy_above_hull_C = calculate_phase_diagram_condition_C(
            all_entries_C,
            entriesGases_C,
            locked_Chem_Potential_C,
            TestMat_entry_C
        )

        # Calculate phase diagram for Condition X
        pd_X, energy_per_atom_X, formation_energy_X, energy_above_hull_X = calculate_phase_diagram_condition_X(
            entriesTotal_X,
            entriesGases_X,
            locked_Chem_Potential_X,
            TestMat_entry_X
        )

        # Output results
        print(f"Results for {comp_str}:")
        print(f"Condition A - Energy Above Hull: {energy_above_hull_A}")
        print(f"Condition C - Energy Above Hull: {energy_above_hull_C}")
        print(f"Condition X - Energy Above Hull: {energy_above_hull_X}")

        # Write results to files
        with open(f'Anode_{comp_str}.txt', 'w') as out_file_A:
            out_file_A.write(f"{comp_str}\n")
            out_file_A.write(f"{energy_above_hull_A}\n")
            out_file_A.write('************************************************************************************************\n')

        with open(f'Cathode_{comp_str}.txt', 'w') as out_file_C:
            out_file_C.write(f"{comp_str}\n")
            out_file_C.write(f"{energy_above_hull_C}\n")
            out_file_C.write('************************************************************************************************\n')

        with open(f'CO2_{comp_str}_X.txt', 'w') as out_file_X:
            out_file_X.write(f"{comp_str}\n")
            out_file_X.write(f"{energy_above_hull_X}\n")
            out_file_X.write('************************************************************************************************\n')


if __name__ == "__main__":
    main()