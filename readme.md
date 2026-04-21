# Protein-Protein Interface (PPI) Predictor

This repository contains an end-to-end pipeline for predicting residue-level protein interfaces. The system integrates structural metrics from **AlphaFold-3**, binding energetics from **PyRosetta**, bioinformatics scores from **BIPSPI**, and evolutionary embeddings from **ESM-2**.

All data and code can be accessed via this link: https://zenodo.org/records/19443030?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImExMTQ4NTk3LWFjYjYtNDA2Mi1hMDk1LTIxNzVmMzlkNmZjYiIsImRhdGEiOnt9LCJyYW5kb20iOiJhYWI4NWRiZWIyNTA5YzVhMTI5MzQwMjJiMDQ5OGJhOCJ9.hLPxJL5nW5Q-OycDz5BVULCjDklY1kKES9dsKQSlNl6R150Y5r63X9QJzUDfs3DPvj4iwu1aJn8-9D5NybDo8g

It is recommended to use the code from the github rather than the Zenodo, as the github will always be the most up-to-date.

## Project Structure

### Required Input Files & Folders
To replicate the results, ensure the following are present in your working directory:

* **`ProteinPairData.py`**: Script used to generate ground-truth labels by downloading and analyzing experimental PDBs.
* **`PPI_Final_Pipeline.ipynb`**: Main notebook for data merging, feature engineering, and model training.
* **`Table_BM5.5.xlsx`**: Protein-Protein Docking Benchmark 5.5 metadata file.
* **`output_binding_energies.csv`**: CSV file containing Rosetta-calculated binding energies (used to extract the `rosetta_dg` column).
* **`/3seed_AF3`**: Directory containing AlphaFold-3 output files. The code deep-scans this folder specifically for `ranked_*.json` and `.zip` files.
* **`/BIPSPI_Results`**: Directory containing the residue-level prediction CSVs from BIPSPI.

### Automatically Generated Files
The scripts will generate the following outputs during execution:

* **`/experimental_pdbs`**: Folder created by `ProteinPairData.py` to store structures downloaded from RCSB.
* **`MASTER_GROUND_TRUTH.csv`**: The output of the labeling script containing the binary interface labels (0 or 1).
* **`DEBUG_REPORT.csv`**: Diagnostic log tracking chain mapping success and residue counts per complex.
* **`FINAL_TRAINING_DATA_WITH_LABELS - Step 3.csv`**: The final integrated dataset merging all biological and structural features.
* **`esm2_embeddings_650M - Step 2.npy`**: The NumPy matrix containing the 1280-dimensional ESM-2 embeddings.
* **`Performance_[TARGET_COMPLEX].png`**: Table visualization showing the ROC-AUC, PR-AUC, and F1-score for the target validation.
* **`Importance_Individual_[TARGET_COMPLEX].png`**: Bar chart visualization of the top 20 individual feature predictors.
* **`Importance_Category_[TARGET_COMPLEX].png`**: Bar chart visualization of feature importance grouped by category (AF3, ESM-2, etc.).

## Execution Workflow

### 1. Label Generation
Run the labeling script first. It reads `Table_BM5.5.xlsx`, downloads the corresponding PDBs into `/experimental_pdbs`, and identifies interface residues via a 5.0Å `NeighborSearch`.

```bash
python ProteinPairData.py
```

### 2. Feature Integration & Training
Open `PPI_Final_Pipeline.ipynb` and run the cells in order:

1.  **Index Results**: The code indexes BIPSPI results, the `output_binding_energies.csv`, and deep-scans the `/3seed_AF3` directory.
2.  **Combine Features**: Aligns pLDDT, iPAE, Rosetta $\Delta G$, and BIPSPI scores with the ground truth labels.
3.  **Generate Embeddings**: Loads the ESM-2 650M model to generate residue embeddings, then reduces them to 32 dimensions via PCA.
4.  **Train Ensemble**: Trains a weighted ensemble of **Random Forest** and **Hist-Gradient Boosting**.
5.  **Error Analysis**: Identifies complexes that were easy or hard to predict. Displays a distribution of PR-AUC results.
6.  **Validate**: Evaluates performance on SARS-CoV-2 targets (**7CZX**, **7D0D**) and generates metrics visualizations.

## Required Libraries
Install the environment using:

```bash
pip install biopython pandas numpy scikit-learn seaborn matplotlib tqdm transformers torch openpyxl requests joblib scipy
```

## Feature Groups
* **AF3 Local/Global**: `af3_plddt`, `af3_ipae_min`, `af3_ipae_mean`, `af3_contact_prob_max`, `iptm`, `ptm`, `ranking_score`.
* **Energetic**: `rosetta_dg` (extracted from `output_binding_energies.csv`).
* **Bioinformatic**: `bipspi_score`.
* **Evolutionary**: ESM-2 650M residue embeddings (PCA-reduced to 32-D).
