from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import sqlite3

DB_PATH = "data/sialic.db"


def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    return {
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "molar_refractivity": Descriptors.MolMR(mol)
    }


def main():
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()

    rows = cur.execute("""
        SELECT compound_id, smiles FROM compounds
    """).fetchall()

    for cid, smiles in rows:

        if not smiles:
            continue

        desc = compute_descriptors(smiles)
        if not desc:
            print(f"Skipped CID {cid}")
            continue

        cur.execute("""
            UPDATE compounds
            SET 
                tpsa = ?,
                rotatable_bonds = ?,
                fsp3 = ?,
                aromatic_rings = ?,
                heavy_atoms = ?,
                formal_charge = ?,
                molar_refractivity = ?
            WHERE compound_id = ?
        """, (
            desc["tpsa"],
            desc["rotatable_bonds"],
            desc["fsp3"],
            desc["aromatic_rings"],
            desc["heavy_atoms"],
            desc["formal_charge"],
            desc["molar_refractivity"],
            cid
        ))

        print(f"Updated CID {cid}")

    con.commit()
    con.close()


if __name__ == "__main__":
    main()