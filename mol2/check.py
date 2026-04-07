import os
from rdkit import Chem

MOL_DIR = r"C:\Users\Harsha\Desktop\Sialic_Acid\mol2"   # change if needed

is_3d_count = 0
is_2d_count = 0

for file in os.listdir(MOL_DIR):
    if not file.endswith(".mol"):
        continue

    path = os.path.join(MOL_DIR, file)

    try:
        mol = Chem.MolFromMolFile(path, removeHs=False)
        if mol is None:
            print(f"❌ Failed: {file}")
            continue

        conf = mol.GetConformer()
        positions = conf.GetPositions()

        # Check if any z != 0
        is_3d = any(abs(pos[2]) > 1e-3 for pos in positions)

        if is_3d:
            print(f"🟢 3D: {file}")
            is_3d_count += 1
        else:
            print(f"🔵 2D: {file}")
            is_2d_count += 1

    except Exception as e:
        print(f"❌ Error {file}: {e}")

print("\n======================")
print(f"3D molecules: {is_3d_count}")
print(f"2D molecules: {is_2d_count}")