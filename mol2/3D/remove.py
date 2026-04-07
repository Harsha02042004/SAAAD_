import os
import re

folder = r"C:\Users\Harsha\Desktop\Sialic_Acid\mol2\3D"

for filename in os.listdir(folder):
    if re.match(r'^\d+\.\s*', filename):
        file_path = os.path.join(folder, filename)
        os.remove(file_path)
        print(f"Deleted: {filename}")