import os
import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors

MOL2_DIR = "mol2/3D"
PNG_DIR = "static/structures/png"
DB = "data/sialic.db"

conn = sqlite3.connect(DB)

success = 0
failed = 0

print("FILES:", os.listdir(MOL2_DIR))

for file in os.listdir(MOL2_DIR):

    if ".mol" not in file.lower():
        continue

    path = os.path.join(MOL2_DIR, file)
    name = os.path.splitext(file)[0]

    png_path = f"structures/png/{name}.png"

    try:
        mol = Chem.MolFromMolFile(path)

        if mol is None:
            print(f"❌ Failed: {file}")
            failed += 1
            continue

        smiles = Chem.MolToSmiles(mol)
        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        logp = Descriptors.MolLogP(mol)

        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

        conn.execute("""
        INSERT OR IGNORE INTO compounds (
            entry_name, smiles,
            molecular_formula, molecular_weight,
            hba, hbd, logp,
            png_path, mol2_path
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            name,
            smiles,
            formula,
            mw,
            hba,
            hbd,
            logp,
            png_path,
            path
        ))

        print(f"✅ {name}")
        success += 1

    except Exception as e:
        print(f"❌ Error {file}: {e}")
        failed += 1

conn.commit()
conn.close()

print("\n====================")
print(f"SUCCESS: {success}")
print(f"FAILED: {failed}")