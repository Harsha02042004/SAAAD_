# regenerate_highlights.py
import os, sqlite3, hashlib, itertools
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

DB_PATH = os.path.join("data", "sialic.db")
OUT_DIR = os.path.join("static", "highlights")
os.makedirs(OUT_DIR, exist_ok=True)

def safe_mol(smi):
    try: return Chem.MolFromSmiles(smi or "")
    except: return None

def draw_highlighted(mol, match, out_path):
    hit_bonds = [b.GetIdx() for b in mol.GetBonds()
                 if b.GetBeginAtomIdx() in match and b.GetEndAtomIdx() in match]
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(mol,
        highlightAtoms=list(match), highlightBonds=hit_bonds,
        highlightAtomColors={i: (0.2, 0.6, 1.0) for i in match},
        highlightBondColors={i: (0.2, 0.6, 1.0) for i in hit_bonds})
    drawer.FinishDrawing()
    with open(out_path, "wb") as f:
        f.write(drawer.GetDrawingText())

con = sqlite3.connect(DB_PATH)
rows = con.execute("SELECT entry_name, smiles FROM compounds").fetchall()
con.close()

compounds = [(name, smi) for name, smi in rows if safe_mol(smi)]
pairs = list(itertools.combinations(compounds, 2))
print(f"Generating highlights for {len(pairs)} pairs...")

for i, ((n1, s1), (n2, s2)) in enumerate(pairs):
    pair    = "_".join(sorted([s1, s2]))
    hash_id = hashlib.md5(pair.encode()).hexdigest()
    path1   = os.path.join(OUT_DIR, f"{hash_id}_1.png")
    path2   = os.path.join(OUT_DIR, f"{hash_id}_2.png")

    if os.path.exists(path1) and os.path.exists(path2):
        continue  # already done

    mol1, mol2 = safe_mol(s1), safe_mol(s2)
    try:
        mcs = rdFMCS.FindMCS([mol1, mol2], timeout=2,
                              atomCompare=rdFMCS.AtomCompare.CompareElements,
                              bondCompare=rdFMCS.BondCompare.CompareOrder)
        if not mcs.smartsString: continue
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None: continue
        draw_highlighted(mol1, mol1.GetSubstructMatch(patt), path1)
        draw_highlighted(mol2, mol2.GetSubstructMatch(patt), path2)
    except Exception as e:
        print(f"  SKIP {n1} vs {n2}: {e}")

    if i % 500 == 0:
        print(f"  {i}/{len(pairs)} done")

print("Done.")