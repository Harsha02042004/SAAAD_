import os
import re
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# =========================
# 📁 PATH SETUP
# =========================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MOL_DIR = os.path.join(BASE_DIR, "mol2", "3D")
PNG_DIR = os.path.join(BASE_DIR, "static", "structures", "png")

os.makedirs(PNG_DIR, exist_ok=True)

# =========================
# 📊 COUNTERS
# =========================
success = 0
failed = 0
skipped = 0

# =========================
# 🧹 SAFE FILENAME
# =========================
def sanitize_filename(name):
    return re.sub(r'[^a-zA-Z0-9_-]', '_', name)

# =========================
# 🧪 LOAD MOLECULE
# =========================
def load_molecule(file_path):
    mol = None

    if file_path.endswith(".mol2"):
        mol = Chem.MolFromMol2File(file_path, sanitize=False)

    if mol is None:
        mol = Chem.MolFromMolFile(file_path, sanitize=False)

    if mol is None:
        return None

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    return mol

# =========================
# 🎨 MAIN LOOP
# =========================
for file in os.listdir(MOL_DIR):

    if not (file.endswith(".mol") or file.endswith(".mol2")):
        continue

    file_path = os.path.join(MOL_DIR, file)
    name = os.path.splitext(file)[0]

    safe_name = sanitize_filename(name)
    png_path = os.path.join(PNG_DIR, f"{safe_name}.png")

    # Skip existing
    if os.path.exists(png_path):
        print(f"⏭️ Skipped: {file}")
        skipped += 1
        continue

    try:
        mol = load_molecule(file_path)

        if mol is None:
            print(f"❌ Failed to load: {file}")
            failed += 1
            continue

        # 🔥 KEY FIX: REMOVE hydrogens (clean look)
        mol = Chem.RemoveHs(mol)

        # 🔥 Generate clean 2D layout
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            print(f"⚠️ 2D coord failed: {file}")

        # 🎨 Draw clean structure
        Draw.MolToFile(
            mol,
            png_path,
            size=(300, 300),
            kekulize=True
        )

        # Verify output
        if not os.path.exists(png_path):
            print(f"❌ PNG not created: {file}")
            failed += 1
            continue

        print(f"✅ Success: {file}")
        success += 1

    except Exception as e:
        print(f"❌ Error {file}: {str(e)}")
        failed += 1

# =========================
# 📈 SUMMARY
# =========================
print("\n======================")
print(f"✅ SUCCESS: {success}")
print(f"❌ FAILED: {failed}")
print(f"⏭️ SKIPPED: {skipped}")
print("======================")