from biopandas.pdb import PandasPdb
import scipy as sp
import pandas as pd
import os
import numpy as np
import Bio
from abnumber import Chain
import math
import enum
import torch
import torch.nn.functional as F
from tqdm import trange
from Bio.PDB import PDBParser
from Bio.PDB.Selection import unfold_entities
from Bio.SeqIO import PdbIO
from Bio.PDB import MMCIFParser, PDBIO
import warnings
import seaborn as sns
from pathlib import Path
from typing import Dict, Tuple, Sequence, List, Optional, Union
from Levenshtein import distance, ratio
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import argparse
from Bio.PDB import MMCIFParser, PDBIO
import tempfile
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from Bio import BiopythonWarning
from pyrosetta import create_score_function
import shutil
import pyrosetta
from Benchmarking.benchmark.ops.all_funcs import *
from Benchmarking.benchmark.ops.protein import *
from Benchmarking.benchmark.ops.benchmark_clean_funcs import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
warnings.simplefilter('ignore', BiopythonWarning)
init('-ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false \
     -no_fconfig -check_cdr_chainbreaks false')
pyrosetta.init()

def cif_to_pose(cif_path, scorefxn=None):
    """Converts a .cif file to a PyRosetta pose."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    
    # Write a temporary PDB for PyRosetta to read
    io = PDBIO()
    io.set_structure(structure)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    io.save(tmp.name)
    
    pose = pose_from_pdb(tmp.name)
    os.unlink(tmp.name)  # clean up temp file
    return pose


def cif_to_temp_pdb(cif_path):
    """Convert a CIF file to a temporary PDB file. Returns the temp PDB path."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, dir=".")
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp.name)
    tmp.close()
    return tmp.name

def ensure_pdb(filepath):
    """If filepath is a .cif, convert to temp PDB and return path. 
    Otherwise return as-is. Also returns a flag for cleanup."""
    if filepath.lower().endswith(".cif"):
        tmp_path = cif_to_temp_pdb(filepath)
        return tmp_path, True  # True = needs cleanup
    return filepath, False
def interface_energy_calc(pred_pdb, pack_scorefxn, prot_type=None):
    chains = list(get_atmseq(pred_pdb))
    print(f"  Chains found: {chains}")

    if len(chains) < 2:
        print("  Skipping: only 1 chain")
        return float('nan')

    # First chain(s) vs last chain(s)
    # For two chains: "A_B"
    # Adjust split if needed for your specific inputs
    chain1 = chains[0]
    chain2 = "".join(chains[1:])
    pred_interface = f"{chain1}_{chain2}"

    print(f"  Interface string: {pred_interface}")

    pred_pose = pose_from_pdb(pred_pdb)
    interface_analyzer = get_interface_analyzer(pred_interface, pack_scorefxn)
    interface_analyzer.apply(pred_pose)
    interface_analyzer_packsep = get_interface_analyzer(pred_interface, pack_scorefxn, pack_sep=True)
    interface_analyzer_packsep.apply(pred_pose)
    binding_energy_dgsep = interface_analyzer_packsep.get_interface_dG()
    return binding_energy_dgsep

def get_interface_analyzer(partner_chain_str, scorefxn, pack_sep=False):
  interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
  interface_analyzer.fresh_instance()
  interface_analyzer.set_pack_input(True)
  interface_analyzer.set_interface(partner_chain_str)
  interface_analyzer.set_scorefunction(scorefxn)
  interface_analyzer.set_compute_interface_energy(True)
  interface_analyzer.set_compute_interface_sc(True)
  interface_analyzer.set_calc_dSASA(True)
  interface_analyzer.set_pack_separated(pack_sep)

  return interface_analyzer

def interface_energy_calc(pred_pdb, pack_scorefxn, prot_type):
    chains = list(get_atmseq(pred_pdb))
    pred_pose = pose_from_pdb(pred_pdb)
    
    if prot_type == 'antibody':
        ab_chains = [x for x in chains if x in ['H', 'L']]
        ag_chains = [x for x in chains if x not in ['H', 'L']]
        if ab_chains and ag_chains:
            pred_interface = f"{''.join(ab_chains)}_{''.join(ag_chains)}"
        else:
            # Fallback: first chain vs rest
            pred_interface = f"{chains[0]}_{''.join(chains[1:])}"
    elif prot_type == 'nanobody':
        nb_chains = [x for x in chains if x in ['H']]
        ag_chains = [x for x in chains if x not in ['H']]
        if nb_chains and ag_chains:
            pred_interface = f"{''.join(nb_chains)}_{''.join(ag_chains)}"
        else:
            pred_interface = f"{chains[0]}_{''.join(chains[1:])}"
    else:
        pred_interface = f"{chains[0]}_{''.join(chains[1:])}"
    
    print(f"  Interface: {pred_interface}")
    interface_analyzer = get_interface_analyzer(pred_interface, pack_scorefxn)
    interface_analyzer.apply(pred_pose)
    interface_analyzer_packsep = get_interface_analyzer(pred_interface, pack_scorefxn, pack_sep=True)
    interface_analyzer_packsep.apply(pred_pose)
    return interface_analyzer_packsep.get_interface_dG()


af3_root = "PATH/OF/AF3FILES"
resultsfilepath = "binding_energies.csv"
default_protein_type = "antibody"

for job_name in sorted(os.listdir(af3_root)):
    job_path = os.path.join(af3_root, job_name)
    if not os.path.isdir(job_path):
        continue

    # Match the naming pattern: fold_{job_name}_model_0.cif
    best_model = f"fold_{job_name}_model_0.cif"
    src = os.path.join(job_path, best_model)

    if os.path.exists(src):
        dst = os.path.join(out_dir, f"{job_name}.cif")
        shutil.copy2(src, dst)
        #print(f"Copied: {src} -> {dst}")
    else:
        print(f"WARNING: {best_model} not found in {job_path}")
        
pdb_dir = out_dir      
# Collect both .pdb and .cif files
struct_files = [f for f in os.listdir(pdb_dir) 
                if f.lower().endswith((".pdb", ".cif"))]

datastructure = pd.DataFrame({
    "AF3_Dir": [os.path.join(pdb_dir, "")] * len(struct_files),
    "AF3_File": struct_files,
    "Protein_type": [default_protein_type] * len(struct_files),
}).reset_index(drop=True)

pdbs = []
prottypes = []
bind_es = []

energy_fxn = "ref2015"
pack_scorefxn = create_score_function(energy_fxn)

for i in trange(datastructure.shape[0]):
    filepath = datastructure.iloc[i].AF3_Dir + datastructure.iloc[i].AF3_File
    prottype = datastructure.iloc[i].Protein_type

    try:
        # Just check chain count for info, then pass the PATH to interface_energy_calc
        pdb_path, is_temp = ensure_pdb(filepath)
        pose = pose_from_pdb(pdb_path)
        print(f"{filepath}: chains={pose.num_chains()}, residues={pose.total_residue()}")

        if pose.num_chains() < 2:
            print(f"  Skipping: only {pose.num_chains()} chain")
            b_E = float('nan')
        else:
            # Pass the temp PDB path, not the pose
            b_E = interface_energy_calc(pdb_path, pack_scorefxn, prottype)

        if is_temp:
            os.unlink(pdb_path)

    except Exception as e:
        print(f"FAILED on {filepath}: {e}")
        b_E = float('nan')

    bind_es.append(b_E)
    pdbs.append(datastructure.iloc[i].AF3_File)
    prottypes.append(prottype)
    
binding_es = pd.DataFrame({
    "Proteins": pdbs,
    "del_G_B": bind_es,
})
binding_es.to_csv(resultsfilepath, index=False)
binding_es.head()


