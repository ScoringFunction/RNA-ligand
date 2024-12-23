import os
import csv
import MDAnalysis as mda
from Bio.PDB import PDBParser
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
import subprocess
import numpy as np
from scipy.spatial.distance import cdist

# Folder path where your PDB files are located
folder_path = 'Complexes'

# Define common ions and molecules to ignore for ligand identification
common_ions = {'NA', 'K', 'MG', 'ZN', 'CA', 'CL', 'MN', 'FE', 'CU', 'CO', 'BR', 'IOD', 'NI', 'HG', 'AG', 'CD', 'AU', 'PB', 'RB', 'HOH'}

# Known SMILES dictionary for frequently missing ligands
known_smiles = {
    "BTN": "O=C1CCC(N2C(C1=O)C(C(=O)NC2=O)C(=O)O)C",
    "ACT": "CC(=O)O",
    "P13": "C1CCCCC1",
    "AP7": "C1=CC=CC=C1",
    "ISI": "CCO",
    "MGR": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "P14": "CC(=O)NCCC(=O)O",
    "HPA": "C1=CC=C(C=C1)O"
}

# Initialize list to store all features
all_features = []

# Initialize PDB parser
parser = PDBParser(QUIET=True)

# Function to extract RNA-specific features
def extract_rna_features(pdb_file, pdb_id):
    features = {}
    structure = parser.get_structure(pdb_id, pdb_file)
    rna_sequence = ""
    valid_residues = {"A", "U", "G", "C"}
    rna_chain_id = None

    for chain in structure.get_chains():
        for residue in chain:
            if residue.resname.strip() in valid_residues:
                rna_chain_id = chain.id
                break
        if rna_chain_id:
            break

    if rna_chain_id:
        for residue in structure[0][rna_chain_id]:
            if residue.resname.strip() in valid_residues:
                rna_sequence += residue.resname.strip()
        
        features['nucleotide_composition_A'] = rna_sequence.count('A') / len(rna_sequence)
        features['nucleotide_composition_U'] = rna_sequence.count('U') / len(rna_sequence)
        features['nucleotide_composition_G'] = rna_sequence.count('G') / len(rna_sequence)
        features['nucleotide_composition_C'] = rna_sequence.count('C') / len(rna_sequence)
        features['gc_content'] = (rna_sequence.count('G') + rna_sequence.count('C')) / len(rna_sequence)

        # Save RNA sequence to a unique file for each PDB
        rna_sequence_file = f"rna_sequence_{pdb_id}.txt"
        with open(rna_sequence_file, "w") as file:
            file.write(rna_sequence)

        # RNAfold MFE calculation (using ViennaRNA)
        try:
            result = subprocess.run(['RNAfold', rna_sequence_file], capture_output=True, text=True)
            output = result.stdout
            mfe_line = output.split('\n')[1]
            mfe = float(mfe_line.split('(')[-1].strip(')'))
            features['minimum_free_energy'] = mfe
        except Exception as e:
            features['minimum_free_energy'] = None

        # Atom count
        atom_count = len(mda.Universe(pdb_file).select_atoms("nucleic"))
        features['atom_count'] = atom_count

    features['pdb_id'] = pdb_id
    return features

# Function to detect ligands in a PDB file
def detect_ligands(pdb_file):
    structure = parser.get_structure("Structure", pdb_file)
    ligands = []
    
    # Loop over all residues in the PDB structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # If the residue is not a standard RNA nucleotide (A, U, G, C) and not an ion
                if (residue.get_resname() not in ['A', 'U', 'G', 'C']) and (residue.get_resname().upper() not in common_ions):
                    ligands.append(residue.get_resname())  # Add the ligand to the list
    
    return ligands

# Function to extract ligand features (SMILES, molecular weight, hydrogen bond donors/acceptors, etc.)
def extract_ligand_features(ligand_name):
    features = {}
    
    try:
        # Try to fetch ligand SMILES from PubChem
        compounds = pcp.get_compounds(ligand_name, 'name')
        if compounds:
            smiles = compounds[0].isomeric_smiles
            print(f"Found SMILES for {ligand_name}: {smiles}")
            
            # Generate RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                features['SMILES'] = smiles
                features['Molecular Weight (g/mol)'] = CalcExactMolWt(mol)
                features['Hydrogen Bond Donors'] = Lipinski.NumHDonors(mol)
                features['Hydrogen Bond Acceptors'] = Lipinski.NumHAcceptors(mol)
                features['LogP'] = Crippen.MolLogP(mol)
                features['TPSA (Å²)'] = rdMolDescriptors.CalcTPSA(mol)
                features['Atom Count'] = mol.GetNumAtoms()
                features['Chiral Atom Count'] = len([atom for atom in mol.GetAtoms() if atom.HasProp("_ChiralityPossible")])
            else:
                print(f"Failed to generate molecule for SMILES: {smiles}")
        else:
            print(f"No compounds found in PubChem for ligand: {ligand_name}")
    
    except Exception as e:
        print(f"Error extracting features for ligand {ligand_name}: {e}")
    
    return features

# Function to extract complex-specific features (electrostatic and van der Waals interactions)
def extract_complex_features(pdb_file, pdb_id):
    features = {}
    u = mda.Universe(pdb_file)
    rna_atoms = u.select_atoms("resname A C G U")
    drug_atoms = u.select_atoms("not resname A C G U")
    
    if len(rna_atoms) == 0 or len(drug_atoms) == 0:
        print(f"No RNA or ligand atoms found in {pdb_id}")  # Debugging line
        return None

    rna_positive_atoms = rna_atoms.select_atoms("name N*")
    drug_negative_atoms = drug_atoms.select_atoms("name O*")

    electrostatic_distances = cdist(rna_positive_atoms.positions, drug_negative_atoms.positions)
    electrostatic_contacts = electrostatic_distances < 5.0
    num_electrostatic_contacts = np.sum(electrostatic_contacts)
    avg_electrostatic_distance = np.mean(electrostatic_distances[electrostatic_contacts]) if np.sum(electrostatic_contacts) > 0 else None

    features['num_electrostatic_contacts'] = num_electrostatic_contacts
    features['avg_electrostatic_distance'] = avg_electrostatic_distance

    distances = cdist(rna_atoms.positions, drug_atoms.positions)
    vdw_contacts = distances < 4.0
    num_vdw_contacts = np.sum(vdw_contacts)
    avg_vdw_distance = np.mean(distances[vdw_contacts]) if np.sum(vdw_contacts) > 0 else None

    features['num_vdw_contacts'] = num_vdw_contacts
    features['avg_vdw_distance'] = avg_vdw_distance

    return features

# Function to ask the user for the processing option
def ask_for_option():
    print("Select an option:")
    print("1: Process all PDB files in the folder")
    print("2: Process a specific PDB file")

    choice = input("Enter 1 or 2: ").strip()

    if choice == "1":
        return None  # No specific PDB, process all
    elif choice == "2":
        specific_pdb = input("Enter the name of the specific PDB file: ").strip()
        return specific_pdb  # Specific PDB file to process
    else:
        print("Invalid option, please choose either 1 or 2.")
        return ask_for_option()

# Ask user for option to run
specific_pdb = ask_for_option()

# Process each PDB file and collect features
for file_name in os.listdir(folder_path):
    if file_name.endswith(".pdb"):
        print(f"\nProcessing PDB file: {file_name}")  # Debugging line
        
        if specific_pdb and file_name != specific_pdb:
            continue  # Skip files if we're processing only a specific one

        pdb_file = os.path.join(folder_path, file_name)
        pdb_id = os.path.splitext(file_name)[0]

        # Step 1: Detect ligands in the PDB file
        ligands = detect_ligands(pdb_file)
        print(f"Detected ligands for {pdb_id}: {ligands}")  # Debugging line

        # Step 2: Extract features for RNA, ligands, and complex interactions
        rna_features = extract_rna_features(pdb_file, pdb_id)
        print(f"RNA Features for {pdb_id}: {rna_features}")  # Debugging line

        ligand_features = {}
        for ligand in ligands:
            ligand_features = extract_ligand_features(ligand)
            print(f"Features for ligand {ligand}: {ligand_features}")  # Debugging line

        complex_features = extract_complex_features(pdb_file, pdb_id)
        print(f"Complex features for {pdb_id}: {complex_features}")  # Debugging line

        # Combine features into a single dictionary and add to list
        if rna_features or ligand_features or complex_features:
            combined_features = {
                'pdb_id': pdb_id,
                **rna_features,
                **ligand_features,
                **complex_features
            }

            # Append to the list of all features
            all_features.append(combined_features)

# Determine the output file name
if all_features:
    if specific_pdb:
        # Save features for specific PDB file with its ID in the filename
        output_filename = f"{specific_pdb}_features.csv"
    else:
        # Save features for all PDB files
        output_filename = "Features.csv"
    
    fieldnames = list(all_features[0].keys())  # Get field names from the first feature set
    with open(output_filename, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_features)

    print(f"\nFeatures have been successfully merged and saved to '{output_filename}'.")
else:
    print("No features extracted.")
