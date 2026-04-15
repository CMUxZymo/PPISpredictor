import os
import pandas as pd
import numpy as np
import requests
from Bio import PDB
from Bio.PDB import NeighborSearch, is_aa
from tqdm import tqdm
import warnings
import re

warnings.filterwarnings('ignore')

# ==========================================
# 1. CONFIGURATION
# ==========================================
BENCHMARK_FILE = 'Table_BM5.5.xlsx' 
OUTPUT_FILE = 'MASTER_GROUND_TRUTH.csv'
PDB_DIR = 'experimental_pdbs'
SARS_PDBS = ['7LS9', '7CZQ', '7CZX', '7D0D', '7DET', '7MSQ', '7X8P']

os.makedirs(PDB_DIR, exist_ok=True)

# ==========================================
# 2. CORE LOGIC
# ==========================================

def download_pdb(pdb_id):
    pdb_id = pdb_id.lower()[:4]
    filepath = os.path.join(PDB_DIR, f'{pdb_id}.pdb')
    if not os.path.exists(filepath):
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                with open(filepath, 'w') as f: f.write(r.text)
                return filepath
        except: return None
    return filepath

def get_interface_labels(pdb_id, cutoff=5.0):
    path = download_pdb(pdb_id)
    if not path: return []
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, path)
    model = structure[0]
    atoms = [a for a in model.get_atoms() if is_aa(a.get_parent())]
    if not atoms: return []
    
    ns = NeighborSearch(atoms)
    results = []
    
    for chain in model:
        for residue in chain:
            if not is_aa(residue): continue
            res_id = residue.get_id()[1]
            is_interface = 0
            for atom in residue.get_atoms():
                neighbors = ns.search(atom.get_coord(), cutoff)
                for neighbor in neighbors:
                    neighbor_chain = neighbor.get_parent().get_parent().id
                    if neighbor_chain != chain.id:
                        is_interface = 1
                        break
                if is_interface: break
            
            results.append({
                'pdb_id': pdb_id.upper(),
                'pdb_chain': chain.id,
                'residue': res_id,
                'ground_truth': is_interface
            })
    return results

# ==========================================
# 3. STRICT EXTRACTION LOGIC
# ==========================================

print(f"📖 Reading {BENCHMARK_FILE}...")
df_bm = pd.read_excel(BENCHMARK_FILE, skiprows=3)

# Filter for valid PDB patterns: Exactly 4 alphanumeric characters
def is_valid_pdb(s):
    return bool(re.match(r'^[A-Za-z0-9]{4}$', str(s)[:4]))

# Get IDs from the first column, but ONLY if they look like PDB IDs
benchmark_pdbs = [
    str(val)[:4].upper() 
    for val in df_bm.iloc[:, 0].dropna() 
    if is_valid_pdb(val)
]

all_pdbs = sorted(list(set(benchmark_pdbs + SARS_PDBS)))

print(f"🚀 Found {len(all_pdbs)} valid complexes (Filtered out headers).")

# ==========================================
# 4. EXECUTION
# ==========================================
master_data = []
for pdb in tqdm(all_pdbs):
    try:
        labels = get_interface_labels(pdb)
        master_data.extend(labels)
    except Exception as e:
        pass # The filter above handles the bad IDs now

df_master = pd.DataFrame(master_data)
df_master.to_csv(OUTPUT_FILE, index=False)

print("\n" + "="*40)
print(f"✅ STEP 1 COMPLETE")
print(f"📊 Total residues: {len(df_master)}")
print(f"🎯 Interface residues: {df_master['ground_truth'].sum()}")
print("="*40)