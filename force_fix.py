from rdkit import Chem
from rdkit.Chem import AllChem

# Use a clean SMILES (approximation of your structure)
smiles = "CC(=O)NC(CO)C(O)C(O)C(O)C(O)CO"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# Generate 3D
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

# Save clean MOL
Chem.MolToMolFile(mol, "mol2/3D/Neu5AcaMe_clean.mol")

print("✅ Clean MOL generated")