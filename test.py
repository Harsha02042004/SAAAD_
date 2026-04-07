import sqlite3
from rdkit import Chem
from rdkit.Chem import MolStandardize

DB_PATH = "data/sialic.db"

def get_largest_fragment(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    try:
        lfc = MolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        return Chem.MolToSmiles(mol)
    except:
        return None


def clean_database():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    rows = cursor.execute("SELECT compound_id, entry_name, smiles FROM compounds").fetchall()

    updated = 0
    removed = 0

    for cid, name, smiles in rows:
        if not smiles:
            continue

        # Detect multi-fragment
        if "." in smiles:
            print(f"[FIXING] {name}")

            cleaned = get_largest_fragment(smiles)

            if cleaned:
                cursor.execute(
                    "UPDATE compounds SET smiles=? WHERE compound_id=?",
                    (cleaned, cid)
                )
                updated += 1
                print(f"  → cleaned: {cleaned}")
            else:
                cursor.execute(
                    "DELETE FROM compounds WHERE compound_id=?",
                    (cid,)
                )
                removed += 1
                print(f"  → removed (invalid)")

    conn.commit()
    conn.close()

    print("\n✅ DONE")
    print(f"Updated: {updated}")
    print(f"Removed: {removed}")


if __name__ == "__main__":
    clean_database()