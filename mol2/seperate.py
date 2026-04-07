import os
import shutil
from rdkit import Chem

BASE_DIR = "."   # current folder
DIR_2D = os.path.join(BASE_DIR, "2D")
DIR_3D = os.path.join(BASE_DIR, "3D")

os.makedirs(DIR_2D, exist_ok=True)
os.makedirs(DIR_3D, exist_ok=True)

count_2d = 0
count_3d = 0
failed = 0

for file in os.listdir(BASE_DIR):
    if not file.endswith(".mol"):
        continue

    path = os.path.join(BASE_DIR, file)

    try:
        mol = Chem.MolFromMolFile(path, removeHs=False)
        if mol is None:
            print(f"❌ Failed: {file}")
            failed += 1
            continue

        conf = mol.GetConformer()
        positions = conf.GetPositions()

        # Check z-coordinates
        is_3d = any(abs(pos[2]) > 1e-3 for pos in positions)

        if is_3d:
            shutil.move(path, os.path.join(DIR_3D, file))
            print(f"🟢 3D → {file}")
            count_3d += 1
        else:
            shutil.move(path, os.path.join(DIR_2D, file))
            print(f"🔵 2D → {file}")
            count_2d += 1

    except Exception as e:
        print(f"❌ Error {file}: {e}")
        failed += 1

print("\n======================")
print(f"2D files: {count_2d}")
print(f"3D files: {count_3d}")
print(f"Failed: {failed}")