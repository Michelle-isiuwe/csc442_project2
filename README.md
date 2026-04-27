# 🧬 DNA & RNA Sequence Analyzer

**CSC 442 — Computational Biology & Interdisciplinary Studies | Project 2**

A production-ready multi-page Streamlit web application for analysing DNA and RNA sequences, performing transcription & translation, searching the UniProt protein database, and visualising sequence statistics.

---

## 📌 Features

### 🔬 Sequence Analyzer
- Paste sequence or upload `.txt` / `.fasta` files (drag-and-drop supported)
- **Auto-detection** of DNA, RNA, or invalid sequences with clear explanation
- **Strand type selection** — Template or Coding strand
- **Transcription** — DNA → mRNA with step-by-step explanation
- **Translation** — mRNA → Codons → Amino Acids (full name, 3-letter, 1-letter)
- **Polypeptide display** with beginner-friendly explanation

### 🧫 Protein Search
- Integrates the **UniProt REST API** for real protein data
- Displays: protein name, organism, function, accession, link to UniProt page
- Results saved to SQLite

### 📊 Visualizations
- Base composition pie chart
- GC vs AT/AU content bar chart
- Amino acid frequency histogram
- Full sequence statistics table

### 📋 History
- All sequence analyses stored in SQLite
- Protein search history
- Search, browse, delete, and export CSV

---

## 🛠️ Installation & Local Setup

```bash
unzip project2_dna_rna_analyzer.zip
cd project2_dna_rna_analyzer

# Optional virtual environment
python3 -m venv venv
source venv/bin/activate      # macOS/Linux
venv\Scripts\activate         # Windows

pip install -r requirements.txt
streamlit run app.py
```

App opens at **http://localhost:8501**.

---

## 🗄️ SQLite Setup

Database (`database/dna_rna.db`) is created automatically on first run. Tables:

- `sequences` — stores input, type, strand, mRNA, protein, stats, timestamp
- `protein_searches` — stores UniProt search queries and results

---

## 🌐 UniProt API

The protein search uses the **UniProt REST API** (`https://rest.uniprot.org`).  
No API key required. An internet connection is needed for this feature.

---

## ☁️ Streamlit Community Cloud Deployment

1. Push the project to a public GitHub repository.
2. Visit [share.streamlit.io](https://share.streamlit.io) and sign in.
3. Click **"New app"**, choose your repo, set main file to `app.py`.
4. Deploy.

> **Note:** SQLite is ephemeral on Streamlit Cloud. For persistence, migrate to a cloud database.

---

## 📁 Folder Structure

```
project2_dna_rna_analyzer/
├── app.py                     # Main multi-page Streamlit app
├── requirements.txt
├── README.md
├── .streamlit/
│   └── config.toml
├── database/
│   ├── __init__.py
│   └── db.py                  # SQLite init, CRUD
├── modules/
│   ├── __init__.py
│   ├── analyzer.py            # Detection, transcription, translation, codon table
│   └── protein_api.py         # UniProt REST API integration
├── utils/
│   ├── __init__.py
│   └── helpers.py             # File saving, CSV export, truncation
├── uploads/                   # Uploaded sequence files
├── assets/
└── tests/
    └── __init__.py
```

---

## 🧪 Supported Sequence Formats

- Plain nucleotide strings: `ATGCGATCG...`
- FASTA format: `>Header\nATGCGATC...`
- `.txt` and `.fasta` file uploads

---

## 📸 Example Screenshots

> *(Add screenshots after first run.)*

---

## 📄 License

For academic use only — CSC 442, Department of Computer Science.
