import pandas as pd
import requests
from Bio import PDB
import numpy as np
import os

def download_pdb(pdb_id, pdb_dir):
    """Download PDB file from RCSB"""
    pdb_id_clean = pdb_id.split('_')[0].lower()
    filepath = os.path.join(pdb_dir, f'{pdb_id_clean}.pdb')
    
    if os.path.exists(filepath):
        return filepath
    
    url = f'https://files.rcsb.org/download/{pdb_id_clean}.pdb'
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            with open(filepath, 'w') as f:
                f.write(response.text)
            return filepath
    except:
        pass
    return None

def extract_sequence(pdb_file, chain_id):
    """Extract amino acid sequence from PDB file for a specific chain"""
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        ppb = PDB.PPBuilder()
        sequences = []
        for pp in ppb.build_peptides(structure[0][chain_id]):
            sequences.append(str(pp.get_sequence()))
        return ''.join(sequences)
    except:
        return None

def extract_interface_residues(bound_pdb, chain1, chain2, threshold=5.0):
    """Compute interface residues (residues within threshold Å of partner)"""
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('complex', bound_pdb)
        chain_a = structure[0][chain1]
        chain_b = structure[0][chain2]
        
        interface_a = set()
        interface_b = set()
        
        for res_a in chain_a.get_residues():
            for atom_a in res_a.get_atoms():
                for res_b in chain_b.get_residues():
                    for atom_b in res_b.get_atoms():
                        if atom_a - atom_b < threshold:
                            interface_a.add(res_a.get_id()[1])
                            interface_b.add(res_b.get_id()[1])
        
        return interface_a, interface_b
    except:
        return set(), set()

def create_labels(sequence, interface_positions):
    """Create binary labels: 1 if residue at interface, 0 otherwise"""
    labels = np.zeros(len(sequence), dtype=int)
    for pos in interface_positions:
        if 1 <= pos <= len(sequence):
            labels[pos - 1] = 1
    return labels

# ============================================================================
# MAIN EXTRACTION
# ============================================================================

benchmark_file = r'C:\Users\aramc\Downloads\Table_BM5.5.xlsx'
num_pairs = 10
output_file = 'benchmark_dataset.csv'
pdb_dir = 'pdb_files'

# Create PDB directory
os.makedirs(pdb_dir, exist_ok=True)

# Read benchmark table
print("Reading benchmark table...")
df = pd.read_excel(benchmark_file, header=0)
df.columns = ['Complex', 'Category', 'PDB_ID_1', 'Protein_1', 'PDB_ID_2', 'Protein_2', 'I_RMSD', 'DASA', 'BM_version']
df_clean = df[df['Complex'].notna() & (df['Complex'].str.contains('_', na=False, regex=True))].copy()

print(f"Found {len(df_clean)} total complexes in benchmark")
print(f"Processing first {min(num_pairs, len(df_clean))} pairs...\n")

# Store results
results = []

# Process N pairs
for idx in range(min(num_pairs, len(df_clean))):
    row = df_clean.iloc[idx]
    complex_id = row['Complex']
    pdb_id_1 = row['PDB_ID_1']
    pdb_id_2 = row['PDB_ID_2']
    bound_id = complex_id.split('_')[0]
    
    # Extract chain IDs from complex notation
    try:
        chain_str = complex_id.split('_')[1]
        if ':' in chain_str:
            chain_1 = chain_str.split(':')[0][0]
            chain_2 = chain_str.split(':')[1][0]
        else:
            chain_1 = chain_str[0]
            chain_2 = chain_str[1] if len(chain_str) > 1 else 'A'
    except:
        chain_1, chain_2 = 'A', 'B'
    
    print(f"[{idx+1}/{min(num_pairs, len(df_clean))}] {complex_id}")
    
    # Download bound complex
    bound_file = download_pdb(bound_id, pdb_dir)
    if bound_file is None:
        print(f"  ✗ Could not download {bound_id}")
        continue
    
    # Extract sequences
    seq1 = extract_sequence(bound_file, chain_1)
    seq2 = extract_sequence(bound_file, chain_2)
    
    if seq1 is None or seq2 is None:
        print(f"  ✗ Could not extract sequences")
        continue
    
    # Extract interface residues
    interface_1, interface_2 = extract_interface_residues(bound_file, chain_1, chain_2)
    
    # Create labels
    labels1 = create_labels(seq1, interface_1)
    labels2 = create_labels(seq2, interface_2)
    
    # Store result
    result = {
        'complex_id': complex_id,
        'protein_a_sequence': seq1,
        'protein_a_length': len(seq1),
        'protein_a_interface_count': int(np.sum(labels1)),
        'protein_a_interface_labels': ','.join(map(str, labels1)),
        'protein_b_sequence': seq2,
        'protein_b_length': len(seq2),
        'protein_b_interface_count': int(np.sum(labels2)),
        'protein_b_interface_labels': ','.join(map(str, labels2))
    }
    
    results.append(result)
    
    print(f"  ✓ Protein A: {len(seq1)} residues, {int(np.sum(labels1))} at interface")
    print(f"  ✓ Protein B: {len(seq2)} residues, {int(np.sum(labels2))} at interface")

# Save to CSV
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print(f"\n✓ Extracted {len(results)} pairs")
print(f"✓ Saved to {output_file}")
print(f"\nDataset summary:")
print(f"  Total residues: {results_df['protein_a_length'].sum() + results_df['protein_b_length'].sum()}")
print(f"  Total interface residues: {results_df['protein_a_interface_count'].sum() + results_df['protein_b_interface_count'].sum()}")