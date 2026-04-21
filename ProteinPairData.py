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
DEBUG_FILE = 'DEBUG_REPORT.csv' # <--- NEW DEBUG OUTPUT
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
    if not path: 
        return [], {"error": "Download failed", "chains": [], "interface_count": 0}
    
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, path)
        model = structure[0]
    except Exception as e:
        return [], {"error": f"Parse error: {str(e)}", "chains": [], "interface_count": 0}

    # Debug: Record the exact chains present in the experimental PDB
    chain_ids = [chain.id for chain in model]
    
    atoms = [a for a in model.get_atoms() if is_aa(a.get_parent())]
    if not atoms: 
        return [], {"error": "No amino acids found", "chains": chain_ids, "interface_count": 0}
    
    ns = NeighborSearch(atoms)
    results = []
    interface_count = 0 # Track how many interfaces we actually find
    
    for chain in model:
        for residue in chain:
            if not is_aa(residue): continue
            res_id = residue.get_id()[1]
            is_interface = 0
            
            for atom in residue.get_atoms():
                neighbors = ns.search(atom.get_coord(), cutoff)
                for neighbor in neighbors:
                    neighbor_chain = neighbor.get_parent().get_parent().id
                    # If an atom is within 5.0A of a DIFFERENT chain, it's an interface
                    if neighbor_chain != chain.id:
                        is_interface = 1
                        break
                if is_interface: break
            
            if is_interface:
                interface_count += 1
                
            results.append({
                'pdb_id': pdb_id.upper(),
                'pdb_chain': chain.id,
                'residue': res_id,
                'ground_truth': is_interface
            })
            
    debug_info = {
        "error": "None",
        "chains": " | ".join(chain_ids),
        "interface_count": interface_count,
        "total_residues": len(results)
    }
    
    return results, debug_info

# ==========================================
# 3. STRICT EXTRACTION LOGIC
# ==========================================

print(f"📖 Reading {BENCHMARK_FILE}...")
try:
    df_bm = pd.read_excel(BENCHMARK_FILE, skiprows=3)
    def is_valid_pdb(s):
        return bool(re.match(r'^[A-Za-z0-9]{4}$', str(s)[:4]))
    benchmark_pdbs = [str(val)[:4].upper() for val in df_bm.iloc[:, 0].dropna() if is_valid_pdb(val)]
except Exception as e:
    print(f"⚠️ Error reading benchmark file: {e}")
    benchmark_pdbs = []

all_pdbs = sorted(list(set(benchmark_pdbs + SARS_PDBS)))
print(f"🚀 Found {len(all_pdbs)} valid complexes to process.")

# ==========================================
# 4. EXECUTION & DEBUG TRACKING
# ==========================================
master_data = []
debug_records = []

for pdb in tqdm(all_pdbs):
    labels, debug_info = get_interface_labels(pdb)
    
    if labels:
        master_data.extend(labels)
        
    debug_records.append({
        'pdb_id': pdb.upper(),
        'chains_in_pdb': debug_info.get('chains', ''),
        'total_residues': debug_info.get('total_residues', 0),
        'interface_residues_found': debug_info.get('interface_count', 0),
        'status': debug_info.get('error', 'Unknown')
    })

# Save the Ground Truth
df_master = pd.DataFrame(master_data)
df_master.to_csv(OUTPUT_FILE, index=False)

# Save the Debug Report
df_debug = pd.DataFrame(debug_records)
df_debug.to_csv(DEBUG_FILE, index=False)

# Identify the problem complexes
zero_interfaces = df_debug[df_debug['interface_residues_found'] == 0]

print("\n" + "="*40)
print(f"✅ STEP 1 COMPLETE")
print(f"📊 Total residues processed: {len(df_master)}")
print(f"🎯 Total Interface residues: {df_master['ground_truth'].sum() if not df_master.empty else 0}")
print("="*40)
print(f"⚠️ DEBUG ALERT: {len(zero_interfaces)} complexes returned 0 interfaces.")
print(f"Check '{DEBUG_FILE}' to see which chains the PDB file used vs. what AlphaFold expected.")
print("="*40)