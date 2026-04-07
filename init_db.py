import sqlite3
import os

DB = "data/sialic.db"

os.makedirs("data", exist_ok=True)


def add_column_if_missing(conn, column_def):
    try:
        conn.execute(f"ALTER TABLE compounds ADD COLUMN {column_def}")
        print(f"➕ Added column: {column_def}")
    except Exception as e:
        if "duplicate column name" in str(e).lower():
            print(f"✔ Exists: {column_def}")
        else:
            print(f"❌ Error adding {column_def}: {e}")


def main():
    conn = sqlite3.connect(DB)

    # 🔹 Create base table (only if not exists)
    conn.execute("""
    CREATE TABLE IF NOT EXISTS compounds (
        compound_id INTEGER PRIMARY KEY AUTOINCREMENT,

        entry_name TEXT UNIQUE,
        smiles TEXT,

        molecular_formula TEXT,
        molecular_weight REAL,

        hba INTEGER,
        hbd INTEGER,
        logp REAL,

        png_path TEXT,
        mol2_path TEXT
    )
    """)

    print("📦 Base table ready")

    # 🔹 Ensure all advanced columns exist
    columns = [
        "tpsa REAL",
        "rotatable_bonds INTEGER",
        "fsp3 REAL",
        "aromatic_rings INTEGER",
        "heavy_atoms INTEGER",
        "formal_charge INTEGER",
        "molar_refractivity REAL"
    ]

    for col in columns:
        add_column_if_missing(conn, col)

    conn.commit()
    conn.close()

    print("\n✅ Database fully initialized & upgraded")


if __name__ == "__main__":
    main()