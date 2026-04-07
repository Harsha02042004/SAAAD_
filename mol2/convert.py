import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolAlign

# =========================
# PATHS
# =========================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_DIR = os.path.join(BASE_DIR, "2D")
OUTPUT_DIR = os.path.join(BASE_DIR, "3D_converted")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("INPUT_DIR:", INPUT_DIR)
print("OUTPUT_DIR:", OUTPUT_DIR)
print("FILES:", os.listdir(INPUT_DIR))

# =========================
# SETTINGS
# =========================
NUM_CONFS = 30
MAX_ITERS = 500
RMSD_THRESHOLD = 0.5

success = 0
failed = 0


# =========================
# CLUSTERING FUNCTION
# =========================
def cluster_conformers(mol, conf_ids):
    dists = []

    for i in range(len(conf_ids)):
        for j in range(i):
            rms = rdMolAlign.GetBestRMS(mol, mol, conf_ids[i], conf_ids[j])
            dists.append(rms)

    clusters = Butina.ClusterData(
        dists,
        len(conf_ids),
        RMSD_THRESHOLD,
        isDistData=True
    )

    return clusters


# =========================
# CHECK 3D
# =========================
def has_3d(mol):
    try:
        conf = mol.GetConformer()
        return conf.Is3D()
    except:
        return False


# =========================
# MAIN LOOP
# =========================
for file in os.listdir(INPUT_DIR):

    if ".mol" not in file.lower():
        continue

    path = os.path.join(INPUT_DIR, file)
    out_path = os.path.join(OUTPUT_DIR, file)

    try:
        mol = Chem.MolFromMolFile(path, sanitize=True)

        if mol is None:
            print(f"❌ Failed to read: {file}")
            failed += 1
            continue

        is_3d = has_3d(mol)

        # =========================
        # CASE 1: EXISTING 3D (ChemSketch etc.)
        # =========================
        if is_3d:
            print(f"🟡 Refining existing 3D: {file}")

            mol = Chem.AddHs(mol)

            try:
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                else:
                    ff = AllChem.UFFGetMoleculeForceField(mol)

                ff.Minimize(maxIts=MAX_ITERS)
                energy = ff.CalcEnergy()

                Chem.MolToMolFile(mol, out_path)

                print(f"🔵 Refined: {file} | Energy: {round(energy,3)}")
                success += 1

            except Exception as e:
                print(f"❌ Optimization failed: {file} | {e}")
                failed += 1

        # =========================
        # CASE 2: PURE 2D → FULL PIPELINE
        # =========================
        else:
            print(f"🟢 Generating 3D: {file}")

            mol = Chem.AddHs(mol)

            # EMBEDDING
            params = AllChem.ETKDGv3()
            params.randomSeed = 42

            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=NUM_CONFS,
                params=params
            )

            if len(conf_ids) == 0:
                print(f"❌ Embedding failed: {file}")
                failed += 1
                continue

            energies = {}

            # OPTIMIZATION
            for cid in conf_ids:
                try:
                    if AllChem.MMFFHasAllMoleculeParams(mol):
                        props = AllChem.MMFFGetMoleculeProperties(mol)
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                    else:
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)

                    ff.Minimize(maxIts=MAX_ITERS)
                    energy = ff.CalcEnergy()

                    energies[cid] = energy

                except Exception:
                    print(f"⚠️ Optimization failed: {file} (conf {cid})")

            if not energies:
                print(f"❌ No valid conformers: {file}")
                failed += 1
                continue

            # CLUSTERING
            clusters = cluster_conformers(mol, list(energies.keys()))

            # SELECT REPRESENTATIVE
            best_candidates = []

            for cluster in clusters:
                cluster_energies = [(cid, energies[cid]) for cid in cluster]
                best_in_cluster = min(cluster_energies, key=lambda x: x[1])
                best_candidates.append(best_in_cluster)

            best_conf, best_energy = min(best_candidates, key=lambda x: x[1])

            # SAVE
            Chem.MolToMolFile(mol, out_path, confId=best_conf)

            print(
                f"🔵 {file} | confs: {len(conf_ids)} | clusters: {len(clusters)} | energy: {round(best_energy,3)}"
            )

            success += 1

    except Exception as e:
        print(f"❌ Error {file}: {e}")
        failed += 1


# =========================
# SUMMARY
# =========================
print("\n======================")
print(f"Converted: {success}")
print(f"Failed: {failed}")