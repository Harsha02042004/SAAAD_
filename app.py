# =========================
# IMPORTS
# =========================
import os
import re
import sqlite3
import tempfile
import subprocess
from collections import defaultdict, Counter
from io import BytesIO

from flask import Flask, render_template, request, jsonify, send_file, after_this_request, Response

from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Scaffolds import MurckoScaffold

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw

# =========================
# APP INIT
# =========================
app = Flask(__name__)
DB_PATH = os.path.join("data", "sialic.db")

gen = GetMorganGenerator(radius=2, fpSize=2048)

CLUSTER_COLORS = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6", "#f39c12"]

PNG_DIR = os.path.join("static", "structures", "png")
os.makedirs(PNG_DIR, exist_ok=True)

# =========================
# BASIC UTILITIES
# =========================
def get_db():
    con = sqlite3.connect(DB_PATH)
    con.row_factory = sqlite3.Row
    return con

def safe_mol(smiles):
    try:
        return Chem.MolFromSmiles(smiles or "")
    except:
        return None

def sanitize_filename(name):
    return re.sub(r'[^a-zA-Z0-9_-]', '_', name or "unknown")

def slugify(name):
    """Decorative URL slug — only used for pretty URLs, never for lookup."""
    return re.sub(r'[^a-zA-Z0-9]+', '-', name or "").strip('-').lower()

def add_png_path(record):
    safe_name = sanitize_filename(record.get("entry_name"))
    record["png_path"] = f"structures/png/{safe_name}.png"
    return record


# =========================
# PNG GENERATION
# Renders a molecule to PNG bytes using RDKit.
# =========================
def render_mol_to_png(smiles, width=300, height=300):
    """Return PNG bytes for the given SMILES, or None on failure."""
    mol = safe_mol(smiles)
    if mol is None:
        return None
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        print(f"[PNG] render failed for SMILES {smiles!r}: {e}")
        return None


# =========================
# DYNAMIC PNG ROUTE
#
# Serves static/structures/png/<safe_name>.png.
# If the file is missing, generates it on-the-fly from the DB
# and saves it to disk so subsequent requests are instant.
# =========================
@app.route("/structures/png/<path:filename>")
def serve_structure_png(filename):
    disk_path = os.path.join("static", "structures", "png", filename)

    # ── 1. Already on disk → serve directly ──────────────────────────────────
    if os.path.exists(disk_path) and os.path.getsize(disk_path) > 100:
        return send_file(disk_path, mimetype="image/png")

    # ── 2. Not on disk → look up by entry_name and render ────────────────────
    # Reverse the sanitize_filename transform is lossy, so we query by
    # comparing sanitize_filename(entry_name) against the requested filename.
    base = filename.replace(".png", "")

    con = get_db()
    rows = con.execute("SELECT entry_name, smiles FROM compounds").fetchall()
    con.close()

    smiles = None
    for row in rows:
        if sanitize_filename(row["entry_name"]) == base:
            smiles = row["smiles"]
            break

    if smiles is None:
        # Return a tiny transparent PNG placeholder so <img> doesn't break
        return _placeholder_png(), 404

    png_bytes = render_mol_to_png(smiles)
    if png_bytes is None:
        return _placeholder_png(), 500

    # Save to disk for next time
    try:
        with open(disk_path, "wb") as f:
            f.write(png_bytes)
    except Exception as e:
        print(f"[PNG] could not cache {disk_path}: {e}")

    return Response(png_bytes, mimetype="image/png")


def _placeholder_png():
    """1×1 transparent PNG — keeps <img> tags from showing broken icons."""
    import base64
    b64 = (
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk"
        "YPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
    )
    return Response(base64.b64decode(b64), mimetype="image/png")


# =========================
# SCAFFOLD
# =========================
def compute_scaffold(smiles):
    mol = safe_mol(smiles)
    if not mol:
        return None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold, isomericSmiles=True)
    except:
        return None


# =========================
# HIGHLIGHT
# =========================
def get_highlight_images(smiles1, smiles2, name1="mol1", name2="mol2"):
    import hashlib

    mol1 = safe_mol(smiles1)
    mol2 = safe_mol(smiles2)
    if mol1 is None or mol2 is None:
        return None, None, 0.0

    try:
        pair    = "_".join(sorted([smiles1, smiles2]))
        hash_id = hashlib.md5(pair.encode()).hexdigest()
        file1   = f"highlights/{hash_id}_1.png"
        file2   = f"highlights/{hash_id}_2.png"
        path1   = os.path.join("static", file1)
        path2   = os.path.join("static", file2)

        os.makedirs(os.path.join("static", "highlights"), exist_ok=True)

        # Serve from disk if already rendered
        if (
            os.path.exists(path1) and os.path.getsize(path1) > 100 and
            os.path.exists(path2) and os.path.getsize(path2) > 100
        ):
            return file1, file2, 0.0

        # Find MCS
        mcs = rdFMCS.FindMCS([mol1, mol2], timeout=2,
                              atomCompare=rdFMCS.AtomCompare.CompareElements,
                              bondCompare=rdFMCS.BondCompare.CompareOrder)
        if not mcs.smartsString:
            return None, None, 0.0

        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            return None, None, 0.0

        match1 = mol1.GetSubstructMatch(patt)
        match2 = mol2.GetSubstructMatch(patt)

        def draw_highlighted(mol, match, out_path):
            hit_bonds = []
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in match:
                    hit_bonds.append(bond.GetIdx())

            drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
            drawer.drawOptions().addStereoAnnotation = True
            drawer.DrawMolecule(
                mol,
                highlightAtoms=list(match),
                highlightBonds=hit_bonds,
                highlightAtomColors={i: (0.2, 0.6, 1.0) for i in match},
                highlightBondColors={i: (0.2, 0.6, 1.0) for i in hit_bonds},
            )
            drawer.FinishDrawing()
            with open(out_path, "wb") as f:
                f.write(drawer.GetDrawingText())

        draw_highlighted(mol1, match1, path1)
        draw_highlighted(mol2, match2, path2)

        overlap = len(match1) / min(mol1.GetNumAtoms(), mol2.GetNumAtoms()) if match1 else 0.0
        return file1, file2, round(overlap, 2)

    except Exception as e:
        print(f"[Highlight] FAILED: {e}")
        return None, None, 0.0
# =========================
# SIMILARITY
# =========================
def compute_similarity_v2(target, compounds, threshold=0.4, top_n=10):
    target_mol = safe_mol(target.get("smiles"))
    if not target_mol:
        return []

    target_fp       = gen.GetFingerprint(target_mol)
    target_scaffold = safe_mol(target.get("scaffold"))

    def extract_props(c):
        return np.array([
            c.get("molecular_weight") or 0,
            c.get("logp") or 0,
            c.get("tpsa") or 0,
            c.get("hbd") or 0,
            c.get("hba") or 0
        ])

    def property_similarity(p1, p2):
        diff  = np.abs(p1 - p2)
        norm  = np.array([500, 5, 200, 10, 10])
        score = np.clip(1 - (diff / norm), 0, 1)
        return float(np.mean(score))

    target_props = extract_props(target)
    results      = []

    for c in compounds:
        if c["compound_id"] == target["compound_id"]:
            continue

        mol = safe_mol(c.get("smiles"))
        if not mol:
            continue

        fp       = gen.GetFingerprint(mol)
        tanimoto = TanimotoSimilarity(target_fp, fp)

        scaffold_sim = 0.0
        scaf2 = safe_mol(c.get("scaffold"))
        if target_scaffold and scaf2:
            try:
                scaffold_sim = TanimotoSimilarity(
                    gen.GetFingerprint(target_scaffold),
                    gen.GetFingerprint(scaf2)
                )
            except:
                pass

        prop_sim  = property_similarity(target_props, extract_props(c))
        mcs_score = 0.0
        try:
            mcs = rdFMCS.FindMCS([target_mol, mol], timeout=1)
            if mcs.smartsString:
                patt    = Chem.MolFromSmarts(mcs.smartsString)
                match1  = target_mol.GetSubstructMatch(patt)
                smaller = min(target_mol.GetNumAtoms(), mol.GetNumAtoms())
                if smaller > 0:
                    mcs_score = len(match1) / smaller
        except:
            pass

        final_score = (
            0.40 * tanimoto +
            0.20 * scaffold_sim +
            0.20 * prop_sim +
            0.20 * mcs_score
        )
        if final_score < threshold:
            continue

        explanation = []
        explanation.append(
            "Highly similar overall structure" if tanimoto > 0.75 else
            "Moderately similar structure"     if tanimoto > 0.5  else
            "Low structural similarity"
        )
        explanation.append(
            "Core scaffold conserved"        if scaffold_sim > 0.7 else
            "Partially similar scaffold"     if scaffold_sim > 0.4 else
            "Different core scaffold"
        )
        explanation.append(
            "Physicochemical properties aligned" if prop_sim > 0.7 else
            "Moderate property similarity"       if prop_sim > 0.5 else
            "Different physicochemical profile"
        )
        explanation.append(
            "Large shared substructure"   if mcs_score > 0.6 else
            "Partial substructure overlap" if mcs_score > 0.3 else
            "Minimal shared substructure"
        )

        img1, img2, overlap = get_highlight_images(
            target.get("smiles"), c.get("smiles"),
            target.get("entry_name"), c.get("entry_name")
        )

        c_copy = dict(c)
        c_copy.update({
            "tanimoto":         round(tanimoto, 3),
            "scaffold_sim":     round(scaffold_sim, 3),
            "prop_sim":         round(prop_sim, 3),
            "mcs_score":        round(mcs_score, 3),
            "final_similarity": round(final_score, 3),
            "explanation":      explanation,
            "highlight_1":      img1,
            "highlight_2":      img2,
            "overlap":          round(mcs_score, 2),
        })
        results.append(c_copy)

    return sorted(results, key=lambda x: x["final_similarity"], reverse=True)[:top_n]


# =========================
# PCA
# =========================
def get_feature_vector(c):
    return [
        c.get("molecular_weight") or 0,
        c.get("logp") or 0,
        c.get("hbd") or 0,
        c.get("hba") or 0,
        c.get("tpsa") or 0,
        c.get("rotatable_bonds") or 0,
        c.get("fsp3") or 0,
    ]

def compute_pca(compounds):
    try:
        X = np.array([get_feature_vector(c) for c in compounds])
        if len(X) < 2:
            return [(0, 0)] * len(compounds)
        X      = StandardScaler().fit_transform(X)
        coords = PCA(n_components=2).fit_transform(X)
        return coords.tolist()
    except:
        return [(0, 0)] * len(compounds)

def compute_clusters(compounds, n_clusters=4):
    try:
        X = np.array([get_feature_vector(c) for c in compounds])
        if len(X) < n_clusters:
            return [0] * len(compounds)
        X      = StandardScaler().fit_transform(X)
        labels = KMeans(n_clusters=n_clusters, random_state=42).fit_predict(X)
        return labels.tolist()
    except:
        return [0] * len(compounds)


# =========================
# SCORING SYSTEM
# =========================
def score_mw(mw):
    if mw is None: return 0
    if 200 <= mw <= 500: return 1
    if 150 <= mw < 200 or 500 < mw <= 650: return 0.7
    return 0.3

def score_logp(logp):
    if logp is None: return 0
    if 0 <= logp <= 3: return 1
    if -1 <= logp < 0 or 3 < logp <= 5: return 0.7
    return 0.3

def score_hbd(hbd):
    if hbd is None: return 0
    if hbd <= 3: return 1
    if hbd <= 5: return 0.8
    if hbd <= 8: return 0.5
    return 0.3

def score_hba(hba):
    if hba is None: return 0
    if hba <= 6:  return 1
    if hba <= 10: return 0.8
    if hba <= 14: return 0.5
    return 0.3

def score_tpsa(tpsa):
    if tpsa is None: return 0
    if tpsa <= 120: return 1
    if tpsa <= 180: return 0.7
    return 0.3

def score_rotatable(rb):
    if rb is None: return 0
    if rb <= 6:  return 1
    if rb <= 10: return 0.7
    return 0.4

def score_fsp3(fsp3):
    if fsp3 is None: return 0
    if fsp3 >= 0.4:  return 1
    if fsp3 >= 0.25: return 0.7
    return 0.3

def compute_base_score(c):
    return (
        0.20 * score_mw(c.get("molecular_weight")) +
        0.20 * score_logp(c.get("logp")) +
        0.15 * score_hbd(c.get("hbd")) +
        0.15 * score_hba(c.get("hba")) +
        0.15 * score_tpsa(c.get("tpsa")) +
        0.10 * score_rotatable(c.get("rotatable_bonds")) +
        0.05 * score_fsp3(c.get("fsp3"))
    )

def compute_penalty(c):
    penalty, reasons = 0.0, []
    if (c.get("molecular_weight") or 0) > 500:
        penalty += 0.1
        reasons.append("High MW")
    return penalty, reasons

def compute_final_score(c):
    base             = compute_base_score(c)
    penalty, reasons = compute_penalty(c)
    final            = round(max(0, min(1, base * (1 - penalty))), 3)
    return final, penalty, reasons


# =========================
# RULES
# =========================
def compute_rules(c):
    mw   = c.get("molecular_weight")
    hba  = c.get("hba")
    hbd  = c.get("hbd")
    logp = c.get("logp")
    if mw is None:
        return {}

    lip = []
    if mw > 500:                           lip.append("MW > 500")
    if logp is not None and logp > 5:      lip.append("LogP > 5")
    if hba and hba > 10:                   lip.append("HBA > 10")
    if hbd and hbd > 5:                    lip.append("HBD > 5")

    return {
        "Lipinski": {
            "status":  "PASS" if len(lip) <= 1 else "FAIL",
            "details": ", ".join(lip) if lip else "OK",
        },
        "Veber": {
            "status":  "PASS" if not (hba and hbd and hba + hbd > 12) else "FAIL",
            "details": "High polarity" if (hba and hbd and hba + hbd > 12) else "OK",
        },
        "Ghose": {
            "status":  "PASS" if 160 <= mw <= 480 else "FAIL",
            "details": "MW outside range" if not (160 <= mw <= 480) else "OK",
        },
        "Egan": {
            "status":  "PASS" if not (logp and logp > 5.88) else "FAIL",
            "details": "High LogP" if (logp and logp > 5.88) else "OK",
        },
        "Muegge": {
            "status":  "PASS" if not (mw > 600) else "FAIL",
            "details": "MW > 600" if mw > 600 else "OK",
        },
    }


# =========================
# CLASSIFICATION
# =========================
def classify_overall(score, rules):
    fails = sum(1 for r in rules.values() if r["status"] == "FAIL")
    if score >= 0.75:
        label = "Drug-like" if fails <= 1 else "Moderate"
    elif score >= 0.5:
        label = "Moderate" if fails <= 2 else "Low"
    else:
        label = "Low"
    return {
        "label":        label,
        "failed_rules": [k for k, r in rules.items() if r["status"] == "FAIL"],
    }


# =========================
# CONFIDENCE
# =========================
def compute_confidence(c, penalty=0):
    fields = [
        c.get("molecular_weight"), c.get("logp"), c.get("hbd"),
        c.get("hba"), c.get("tpsa"), c.get("rotatable_bonds"), c.get("fsp3"),
    ]
    available = sum(1 for f in fields if f is not None)
    return round((available / len(fields)) * (1 - penalty), 2)

def confidence_label(conf):
    if conf >= 0.75: return "High"
    elif conf >= 0.5: return "Medium"
    return "Low"


# =========================
# RADAR NORMALIZATION
# =========================
def normalize_for_radar(c):
    return {
        "MW":   min((c.get("molecular_weight") or 0) / 600, 1),
        "LogP": min(abs(c.get("logp") or 0) / 5, 1),
        "HBD":  min((c.get("hbd") or 0) / 10, 1),
        "HBA":  min((c.get("hba") or 0) / 15, 1),
        "TPSA": min((c.get("tpsa") or 0) / 200, 1),
        "Flex": min((c.get("rotatable_bonds") or 0) / 15, 1),
        "3D":   c.get("fsp3") or 0,
    }


# =========================
# INTELLIGENCE LAYER
# =========================
def classify_solubility(c):
    logp = c.get("logp") or 0
    tpsa = c.get("tpsa") or 0
    if logp < 2 and tpsa > 60: return "Good"
    elif logp < 4:              return "Moderate"
    return "Poor"

def classify_permeability(c):
    tpsa = c.get("tpsa") or 0
    if tpsa < 90:  return "High"
    elif tpsa < 140: return "Moderate"
    return "Low"

def classify_bbb(c):
    logp = c.get("logp") or 0
    tpsa = c.get("tpsa") or 0
    return "Likely permeant" if tpsa < 90 and logp > 1 else "Unlikely"

def classify_flexibility(c):
    rb = c.get("rotatable_bonds") or 0
    if rb <= 5:  return "Rigid"
    elif rb <= 10: return "Moderate"
    return "Flexible"

def classify_size(c):
    mw = c.get("molecular_weight") or 0
    if mw < 250:  return "Fragment-like"
    elif mw <= 500: return "Drug-like"
    return "Large"

def generate_assessment(c):
    positives, concerns = [], []
    logp = c.get("logp") or 0
    tpsa = c.get("tpsa") or 0
    mw   = c.get("molecular_weight") or 0
    rb   = c.get("rotatable_bonds") or 0

    if tpsa < 90:      positives.append("High permeability")
    if 0 <= logp <= 3: positives.append("Balanced lipophilicity")
    if mw <= 500:      positives.append("Drug-like molecular size")
    if rb <= 10:       positives.append("Acceptable flexibility")
    if tpsa >= 90:     concerns.append("Limited BBB penetration")
    if logp > 4:       concerns.append("High lipophilicity")
    if rb > 10:        concerns.append("High flexibility")
    return positives, concerns

def generate_summary(c):
    score = c.get("score", 0)
    if score > 0.8:   return "Strong drug-like candidate with favorable properties"
    elif score > 0.6: return "Moderately drug-like compound with some limitations"
    return "Low drug-likeness; requires optimization"

def explain_score(c):
    reasons = []
    if (c.get("molecular_weight") or 0) <= 500: reasons.append("Favorable molecular weight")
    if (c.get("logp") or 0) <= 5:               reasons.append("Balanced lipophilicity")
    if (c.get("tpsa") or 0) <= 140:             reasons.append("Good permeability")
    if (c.get("hbd") or 0) <= 5:                reasons.append("Acceptable HBD")
    if (c.get("hba") or 0) <= 10:               reasons.append("Acceptable HBA")
    return reasons


# =========================
# DOWNLOAD ROUTE
# =========================
THREED_DIR = r"C:\Users\Harsha\Desktop\Sialic_Acid\mol2\3D"

@app.route("/download/<fmt>/<int:cid>")
def download_structure(fmt, cid):
    con = get_db()
    row = con.execute("SELECT * FROM compounds WHERE compound_id=?", (cid,)).fetchone()
    con.close()

    if not row:
        return "Not found", 404

    raw_name  = row["entry_name"]
    safe_name = sanitize_filename(raw_name)

    mol_path = None
    for file in os.listdir(THREED_DIR):
        if file.lower().startswith(raw_name.lower()):
            mol_path = os.path.join(THREED_DIR, file)
            break

    if not mol_path:
        return f"3D file not found for {raw_name}", 404

    tmp_dir = tempfile.mkdtemp()
    try:
        filepath = os.path.join(tmp_dir, f"{safe_name}.{fmt}")
        subprocess.run(["obabel", mol_path, "-O", filepath], check=True)
    except Exception as e:
        return f"Conversion failed: {str(e)}", 500

    @after_this_request
    def cleanup(response):
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return response

    return send_file(filepath, as_attachment=False)


# =========================
# CLUSTER SUMMARY
# =========================
def interpret_cluster(summary):
    if summary["MW"] > 450:                         return "Bulky glycoconjugates"
    elif summary["TPSA"] > 180:                     return "Highly polar surface-binding analogs"
    elif summary["MW"] < 350 and summary["RB"] < 6: return "Compact sialic scaffolds"
    return "Intermediate derivatives"

def summarize_cluster(compounds, cluster_id):
    subset = [c for c in compounds if c.get("cluster") == cluster_id]
    if not subset:
        return {}

    def avg(key):
        vals = [c.get(key) for c in subset if c.get(key) is not None]
        return round(sum(vals) / len(vals), 2) if vals else 0

    summary = {
        "MW":   avg("molecular_weight"),
        "LogP": avg("logp"),
        "TPSA": avg("tpsa"),
        "RB":   avg("rotatable_bonds"),
    }
    summary["label"] = interpret_cluster(summary)
    summary["color"] = CLUSTER_COLORS[cluster_id % len(CLUSTER_COLORS)]
    return summary

def group_pca_points(compounds, precision=3):
    grouped = defaultdict(list)
    for c in compounds:
        x = round(c.get("pca_x", 0), precision)
        y = round(c.get("pca_y", 0), precision)
        grouped[(x, y)].append({
            "compound_id": c["compound_id"],
            "slug":        c.get("slug", slugify(c["entry_name"])),
            "name":        c["entry_name"],
            "score":       c.get("score", 0),
            "cluster":     c.get("cluster", 0),
        })

    result = []
    for (x, y), members in grouped.items():
        cluster_counts = Counter(m["cluster"] for m in members)
        cluster_id     = cluster_counts.most_common(1)[0][0]
        result.append({
            "x":         x,
            "y":         y,
            "count":     len(members),
            "compounds": members,
            "cluster":   cluster_id,
            "color":     CLUSTER_COLORS[cluster_id % len(CLUSTER_COLORS)],
        })
    return result


# =========================
# HELPER: build compound list from DB rows
# =========================
def build_compounds(rows):
    compounds = []
    for r in rows:
        r             = dict(r)
        r             = add_png_path(r)
        r["scaffold"] = compute_scaffold(r.get("smiles"))
        r["slug"]     = slugify(r.get("entry_name", ""))
        score, _, _   = compute_final_score(r)
        r["score"]    = score
        compounds.append(r)
    return compounds


# =========================
# HELPER: paginate
# =========================
def paginate(items, page_str, per_page_str):
    total = len(items)

    if per_page_str == "all":
        return items, 1, 1, "all"

    per_page_int = int(per_page_str)
    total_pages  = max(1, (total + per_page_int - 1) // per_page_int)
    page         = max(1, min(int(page_str), total_pages))
    start        = (page - 1) * per_page_int

    return items[start:start + per_page_int], page, total_pages, per_page_str


# =========================
# HELPER: sort compounds
# =========================
def sort_compounds(compounds, sort):
    if sort == "score_desc":  compounds.sort(key=lambda x: x.get("score", 0), reverse=True)
    elif sort == "score_asc": compounds.sort(key=lambda x: x.get("score", 0))
    elif sort == "mw_desc":   compounds.sort(key=lambda x: x.get("molecular_weight") or 0, reverse=True)
    elif sort == "mw_asc":    compounds.sort(key=lambda x: x.get("molecular_weight") or 0)
    elif sort == "name_asc":  compounds.sort(key=lambda x: (x.get("entry_name") or "").lower())
    elif sort == "name_desc": compounds.sort(key=lambda x: (x.get("entry_name") or "").lower(), reverse=True)
    return compounds


# =========================
# HOMEPAGE
# =========================
@app.route("/")
def homepage():
    con   = get_db()
    total = con.execute("SELECT COUNT(*) FROM compounds").fetchone()[0]
    con.close()
    return render_template("homepage.html", total=total)


# =========================
# SEARCH
# =========================
@app.route("/search")
def search():
    conn = get_db()

    query    = request.args.get("q", "").strip()
    smiles   = request.args.get("smiles", "").strip()
    mw_min   = request.args.get("mw_min")
    mw_max   = request.args.get("mw_max")
    logp_min = request.args.get("logp_min")
    logp_max = request.args.get("logp_max")
    hbd_max  = request.args.get("hbd_max")
    hba_max  = request.args.get("hba_max")
    tpsa_max = request.args.get("tpsa_max")

    similarity_threshold = float(request.args.get("sim", 0.5))
    sort         = request.args.get("sort", "score_desc")
    page_str     = request.args.get("page", "1")
    per_page_str = request.args.get("per_page", "all")

    conditions, params = [], []
    if query:
        conditions.append("(entry_name LIKE ? OR molecular_formula LIKE ?)")
        params.extend([f"%{query}%", f"%{query}%"])
    if mw_min and mw_max:
        conditions.append("molecular_weight BETWEEN ? AND ?")
        params.extend([mw_min, mw_max])
    if logp_min and logp_max:
        conditions.append("logp BETWEEN ? AND ?")
        params.extend([logp_min, logp_max])
    if hbd_max:
        conditions.append("hbd <= ?"); params.append(hbd_max)
    if hba_max:
        conditions.append("hba <= ?"); params.append(hba_max)
    if tpsa_max:
        conditions.append("tpsa <= ?"); params.append(tpsa_max)

    sql = "SELECT * FROM compounds"
    if conditions:
        sql += " WHERE " + " AND ".join(conditions)

    rows      = conn.execute(sql, params).fetchall()
    compounds = build_compounds(rows)

    if smiles:
        query_mol = Chem.MolFromSmiles(smiles)
        if query_mol:
            fpgen    = GetMorganGenerator(radius=2, fpSize=2048)
            query_fp = fpgen.GetFingerprint(query_mol)
            results  = []
            for c in compounds:
                mol = safe_mol(c.get("smiles"))
                if not mol:
                    continue
                sim = TanimotoSimilarity(query_fp, fpgen.GetFingerprint(mol))
                if sim >= similarity_threshold:
                    c["similarity"] = round(sim, 3)
                    results.append(c)
        else:
            results = compounds
    else:
        results = compounds

    if smiles and sort == "score_desc":
        results.sort(key=lambda x: x.get("similarity", 0), reverse=True)
    else:
        sort_compounds(results, sort)

    total = len(results)
    paginated_results, page, total_pages, per_page_str = paginate(results, page_str, per_page_str)
    conn.close()

    return render_template(
        "search.html",
        compounds=paginated_results,
        total=total, page=page, total_pages=total_pages, per_page=per_page_str,
        sort=sort, query=query, smiles=smiles, similarity_threshold=similarity_threshold,
        mw_min=mw_min, mw_max=mw_max, logp_min=logp_min, logp_max=logp_max,
        hbd_max=hbd_max, hba_max=hba_max, tpsa_max=tpsa_max,
    )


# =========================
# API FILTER
# =========================
@app.route("/api/filter")
def api_filter():
    con = get_db()

    query    = request.args.get("q", "").strip()
    mw_min   = request.args.get("mw_min")
    mw_max   = request.args.get("mw_max")
    logp_min = request.args.get("logp_min")
    logp_max = request.args.get("logp_max")
    hbd_max  = request.args.get("hbd_max")
    hba_max  = request.args.get("hba_max")
    tpsa_max = request.args.get("tpsa_max")

    conditions, params = [], []
    if query:
        conditions.append("(entry_name LIKE ? OR molecular_formula LIKE ?)")
        params.extend([f"%{query}%", f"%{query}%"])
    if mw_min and mw_max:
        conditions.append("molecular_weight BETWEEN ? AND ?"); params.extend([mw_min, mw_max])
    if logp_min and logp_max:
        conditions.append("logp BETWEEN ? AND ?"); params.extend([logp_min, logp_max])
    if hbd_max:  conditions.append("hbd <= ?"); params.append(hbd_max)
    if hba_max:  conditions.append("hba <= ?"); params.append(hba_max)
    if tpsa_max: conditions.append("tpsa <= ?"); params.append(tpsa_max)

    sql = "SELECT * FROM compounds"
    if conditions:
        sql += " WHERE " + " AND ".join(conditions)

    rows = con.execute(sql, params).fetchall()
    con.close()
    return jsonify(build_compounds(rows)[:100])


# =========================
# BROWSE
# =========================
@app.route("/browse")
def browse():
    conn = get_db()

    sort         = request.args.get("sort", "score_desc")
    page_str     = request.args.get("page", "1")
    per_page_str = request.args.get("per_page", "all")

    rows      = conn.execute("SELECT * FROM compounds").fetchall()
    compounds = build_compounds(rows)
    sort_compounds(compounds, sort)

    total = len(compounds)
    paginated, page, total_pages, per_page_str = paginate(compounds, page_str, per_page_str)
    conn.close()

    return render_template(
        "browse.html",
        compounds=paginated, total=total, page=page,
        total_pages=total_pages, per_page=per_page_str, sort=sort,
    )


# =========================
# COMPOUND DETAIL
# URL: /compound/<int:cid>/<path:slug>
# cid  → DB lookup  |  slug → decorative only
# =========================
@app.route("/compound/<int:cid>/<path:slug>")
def compound_detail(cid, slug):
    con       = get_db()
    threshold = float(request.args.get("threshold", 0.4))

    rows = con.execute("""
        SELECT compound_id, entry_name, smiles, molecular_formula,
               molecular_weight, logp, hbd, hba, tpsa, rotatable_bonds, fsp3
        FROM compounds
    """).fetchall()
    con.close()

    compounds = build_compounds(rows)

    # ── PCA ──────────────────────────────────────────────────────────────────
    coords = compute_pca(compounds)
    for i, comp in enumerate(compounds):
        comp["pca_x"] = coords[i][0]
        comp["pca_y"] = coords[i][1]

    # ── Cluster + remap left→right by PCA x ──────────────────────────────────
    raw_clusters = compute_clusters(compounds)
    for i, comp in enumerate(compounds):
        comp["cluster"] = raw_clusters[i]

    cluster_centers  = {}
    for cid_k in set(raw_clusters):
        subset = [c for c in compounds if c["cluster"] == cid_k]
        if subset:
            cluster_centers[cid_k] = sum(c["pca_x"] for c in subset) / len(subset)

    ordered_clusters = sorted(cluster_centers, key=lambda k: cluster_centers[k])
    remap            = {old: new for new, old in enumerate(ordered_clusters)}
    for comp in compounds:
        comp["cluster"] = remap[comp["cluster"]]

    # ── PCA bounds ────────────────────────────────────────────────────────────
    xs = [c["pca_x"] for c in compounds]
    ys = [c["pca_y"] for c in compounds]
    pca_bounds = {"x_min": min(xs), "x_max": max(xs), "y_min": min(ys), "y_max": max(ys)}

    # ── Rank ──────────────────────────────────────────────────────────────────
    compounds = sorted(compounds, key=lambda x: x["score"], reverse=True)
    for i, comp in enumerate(compounds):
        comp["rank"] = i + 1

    total_compounds = len(compounds)

    # ── Lookup target ─────────────────────────────────────────────────────────
    target = next((x for x in compounds if x["compound_id"] == cid), None)
    if not target:
        return render_template("404.html"), 404

    # ── PCA groups + cluster summaries ────────────────────────────────────────
    pca_grouped     = group_pca_points(compounds)
    cluster_summary = {i: summarize_cluster(compounds, i) for i in range(len(ordered_clusters))}

    # ── Similarity ────────────────────────────────────────────────────────────
    similar_compounds = compute_similarity_v2(target, compounds, threshold=threshold)

    # ── Scoring + analysis ────────────────────────────────────────────────────
    score, penalty, reasons = compute_final_score(target)
    rules        = compute_rules(target)
    result       = classify_overall(score, rules)
    confidence   = compute_confidence(target, penalty)
    conf_label   = confidence_label(confidence)
    radar        = normalize_for_radar(target)
    explain      = explain_score(target)
    solubility   = classify_solubility(target)
    permeability = classify_permeability(target)
    bbb          = classify_bbb(target)
    flexibility  = classify_flexibility(target)
    size         = classify_size(target)
    positives, concerns = generate_assessment(target)
    summary      = generate_summary(target)

    def get_range(key):
        vals = [c.get(key) for c in compounds if c.get(key) is not None]
        return (min(vals), max(vals)) if vals else (0, 1)

    ranges = {
        "mw":   get_range("molecular_weight"),
        "logp": get_range("logp"),
        "hbd":  get_range("hbd"),
        "hba":  get_range("hba"),
        "tpsa": get_range("tpsa"),
        "rb":   get_range("rotatable_bonds"),
        "fsp3": get_range("fsp3"),
    }

    top_compounds = compounds[:max(5, int(0.1 * len(compounds)))]

    def avg(key):
        vals = [c.get(key) for c in top_compounds if c.get(key) is not None]
        return sum(vals) / len(vals) if vals else 0

    ref_raw = {k: avg(db_key) for k, db_key in [
        ("mw","molecular_weight"),("logp","logp"),("hbd","hbd"),
        ("hba","hba"),("tpsa","tpsa"),("rb","rotatable_bonds"),("fsp3","fsp3"),
    ]}

    def normalize(val, minv, maxv):
        return 0 if maxv == minv else (val - minv) / (maxv - minv)

    reference_profile = [
        normalize(ref_raw["mw"],   *ranges["mw"]),
        normalize(ref_raw["logp"], *ranges["logp"]),
        normalize(ref_raw["hbd"],  *ranges["hbd"]),
        normalize(ref_raw["hba"],  *ranges["hba"]),
        normalize(ref_raw["tpsa"], *ranges["tpsa"]),
        normalize(ref_raw["rb"],   *ranges["rb"]),
        normalize(ref_raw["fsp3"], *ranges["fsp3"]),
    ]

    return render_template(
        "compound.html",
        c=target, rules=rules, score=score, result=result,
        penalty_reasons=reasons, confidence=confidence, confidence_label=conf_label,
        radar=radar, similar_compounds=similar_compounds, threshold=threshold,
        all_compounds=compounds, pca_grouped=pca_grouped, explain=explain,
        solubility=solubility, permeability=permeability, bbb=bbb,
        flexibility=flexibility, size=size, positives=positives, concerns=concerns,
        summary=summary, cluster_summary=cluster_summary, total_compounds=total_compounds,
        CLUSTER_COLORS=CLUSTER_COLORS, pca_bounds=pca_bounds, ranges=ranges,
        reference_profile=reference_profile,
    )


# =========================
# FAQ / CONTACT
# =========================
@app.route("/faq")
def faq():
    return render_template("faq.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")


# =========================
# RUN
# =========================
if __name__ == "__main__":
    app.run(debug=True)