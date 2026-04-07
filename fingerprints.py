from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import sqlite3

DB_PATH = r"C:\Users\Harsha\Desktop\Sialic_Acid\data\sialic.db"

BATCH_SIZE = 1000

# 👉 New RDKit generator (replaces deprecated function)
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def generate_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Invalid SMILES: {smiles}")
        return None

    fp = fpgen.GetFingerprint(mol)
    return fp.ToBitString()


def main():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT compound_id, smiles FROM compounds WHERE fingerprint IS NULL"
    )
    rows = cursor.fetchall()

    print(f"🔍 Found {len(rows)} compounds to process")

    updates = []

    for i, (cid, smiles) in enumerate(rows, start=1):
        fp = generate_fingerprint(smiles)

        if fp:
            updates.append((fp, cid))

        if len(updates) >= BATCH_SIZE:
            cursor.executemany(
                "UPDATE compounds SET fingerprint=? WHERE compound_id=?",
                updates
            )
            conn.commit()
            print(f"✅ Processed {i} rows")
            updates = []

    if updates:
        cursor.executemany(
            "UPDATE compounds SET fingerprint=? WHERE compound_id=?",
            updates
        )
        conn.commit()

    conn.close()
    print("🎉 Fingerprints generation complete!")


if __name__ == "__main__":
    main()