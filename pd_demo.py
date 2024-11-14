from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

# Your material formula
TestMat_formula_A = "Ba8Fe8O24"

# Fetch entries for Ba-Fe-O system and your material
with MPRester("kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK") as m:
    entries = m.get_entries_in_chemsys(["Ba", "Fe", "O"])
    TestMat_entries_A = m.get_entries(TestMat_formula_A)
    if not TestMat_entries_A:
        print(f"No entries found for {TestMat_formula_A}")
        exit()
    TestMat_entry_A = TestMat_entries_A[0]  # Use the first entry found

# Create the phase diagram
pd_A = PhaseDiagram(entries)

# Calculate energy per atom
energy_per_atom = TestMat_entry_A.energy / TestMat_entry_A.composition.num_atoms

# Get hull energy per atom
hull_energy = pd_A.get_hull_energy_per_atom(TestMat_entry_A.composition)
print(f"Hull energy: {hull_energy} eV/atom")

# Compute energy above hull
energy_above_hull_A = energy_per_atom - hull_energy

print(f"Energy above hull: {energy_above_hull_A} eV/atom")