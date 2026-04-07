# 🧬 SAAAAD — Sialic Acid Analogs ADME Prediction Database

A comprehensive, interactive **cheminformatics platform** for exploring, analyzing, and visualizing **sialic acid derivatives** using molecular descriptors, structural similarity, and drug-likeness scoring.

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![Flask](https://img.shields.io/badge/Flask-WebApp-green.svg)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-orange.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Active-success.svg)

---

## 🚀 Overview

SAAAAD (Sialic Acid Analogs ADME Prediction Database) is a Flask-based web application designed to:

* Store and manage **sialic acid derivative compounds**
* Perform **descriptor-based analysis**
* Enable **similarity search (Tanimoto)**
* Visualize **chemical space (PCA)**
* Evaluate **drug-likeness rules**
* Provide **interactive molecular exploration**

This project integrates **RDKit-based cheminformatics pipelines** with a modern web interface for research and discovery.

---

## 🧠 Key Features

### 🔬 Chemical Data & Structures

* SMILES-based compound storage
* 2D molecular structure visualization
* On-demand structure generation (RDKit)

### 📊 Molecular Descriptors

* Molecular weight, LogP, H-bond donors/acceptors
* Topological polar surface area (TPSA)
* Rotatable bonds

### 💊 Drug-Likeness Scoring

* Lipinski Rule of Five
* Veber Rule
* Ghose / Egan / Muegge (partial implementations)
* Custom scoring system

### 🔗 Similarity Analysis

* Morgan fingerprints
* Tanimoto similarity
* Nearest-neighbor identification

### 📈 Chemical Space Visualization

* PCA-based dimensionality reduction
* Cluster visualization
* Compound highlighting

### 🔍 Interactive UI

* Live search filtering
* Pagination & sorting
* Compound detail pages
* Descriptor dashboards

---

## 🏗️ Project Structure

```
SAAAD/
│
├── app.py                  # Main Flask application
├── templates/             # HTML frontend
├── static/                # CSS, JS, images
├── data/                  # Raw dataset
│
├── scripts/               # Preprocessing pipeline
│   ├── init_db.py
│   ├── compute_descriptors.py
│   ├── fingerprints.py
│   └── clean_smiles.py
│
├── requirements.txt
└── README.md
```

---

## ⚙️ Installation

### 1. Clone the repository

```
git clone https://github.com/Harsha02042004/SAAAD_.git
cd SAAAD_
```

### 2. Create environment (recommended)

```
python -m venv venv
venv\Scripts\activate   # Windows
```

### 3. Install dependencies

```
pip install -r requirements.txt
```

---

## ▶️ Running the App

```
python app.py
```

Then open:

```
http://127.0.0.1:5000
```

---

## 🧪 Data Processing Pipeline

The database is built using modular preprocessing scripts:

* `init_db.py` → Initializes database schema
* `clean_smiles.py` → Removes salts & invalid fragments
* `compute_descriptors.py` → Calculates molecular descriptors
* `fingerprints.py` → Generates fingerprints for similarity

⚠️ Large assets (e.g., generated images) are excluded from GitHub and created dynamically or offline.

---

## 📦 Technologies Used

* **Python (Flask)** — Backend framework
* **RDKit** — Cheminformatics engine
* **SQLite** — Lightweight database
* **Chart.js** — Data visualization
* **3Dmol.js** — Molecular visualization

---

## 🎯 Applications

* Drug discovery (glycan-targeted therapeutics)
* Structure–activity relationship (SAR) studies
* Chemical space exploration
* Bioinformatics + cheminformatics integration

---

## 🚧 Future Improvements

* Advanced clustering (scaffold-based)
* ML-based activity prediction
* API endpoints for programmatic access
* Enhanced 3D visualization
* Cloud deployment (Render / Docker)

---

## ⚠️ System Dependencies

This project requires Open Babel for file conversion:

Install Open Babel before running:
- https://openbabel.org/wiki/Main_Page

## ⭐ Contribution & Support

If you find this useful:

* ⭐ Star the repo
* 🍴 Fork and extend
* 🧠 Suggest features

---

## 📜 License

This project is open-source and available under the MIT License.
