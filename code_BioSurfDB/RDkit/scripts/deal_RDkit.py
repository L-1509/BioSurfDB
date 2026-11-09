"""
Combine SMILES from two possible columns (SMILES_PubChem first, then SMILES_Original).
If both are empty, leave chemical properties as NaN and count such rows.
Outputs standardized chemical information using RDKit.
"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import numpy as np

def process_smiles_data(input_file, output_file):
    """
    Process Excel data to generate standardized chemical info from SMILES.

    Priority: SMILES_PubChem > SMILES_Original (if the former is empty or missing).
    Adds columns: SMILES_RDKit, SMILES_Parent_RDKit, InChI_RDKit, InChIKey_RDKit,
                  Molecular_Formula_RDKit, Molecular_Weight_RDKit.

    Parameters:
    -----------
    input_file : str
        Path to input Excel file.
    output_file : str
        Path to output Excel file.
    """
    
    # Check input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Read Excel, first row as header (sheet_name=0 reads the first sheet by position)
    df = pd.read_excel(input_file, sheet_name=0, header=0)
    
    # Detect columns
    col_pubchem = None
    col_original = None
    for col in df.columns:
        if 'SMILES_PubChem' in str(col):
            col_pubchem = col
        if 'SMILES_Original' in str(col):
            col_original = col
    
    if col_pubchem is None and col_original is None:
        raise ValueError("Neither 'SMILES_PubChem' nor 'SMILES_Original' column found.")
    
    print(f"Primary column: {col_pubchem if col_pubchem else 'None'}")
    print(f"Fallback column: {col_original if col_original else 'None'}")
    
    # Storage lists
    clean_smiles_list = []
    parent_smiles_list = []
    inchi_list = []
    inchikey_list = []
    formula_list = []
    mw_list = []
    
    # Counter for rows where both SMILES columns are empty
    both_empty_count = 0
    
    # Process each row
    for idx, row in df.iterrows():
        smiles = None
        # Try primary column first
        if col_pubchem is not None:
            val = row[col_pubchem]
            if pd.notna(val) and str(val).strip() != '':
                smiles = str(val).strip()
        # If still None, try fallback
        if smiles is None and col_original is not None:
            val = row[col_original]
            if pd.notna(val) and str(val).strip() != '':
                smiles = str(val).strip()
        
        # If both are empty, count and set NaN
        if smiles is None:
            both_empty_count += 1
            clean_smiles_list.append(np.nan)
            parent_smiles_list.append(np.nan)
            inchi_list.append(np.nan)
            inchikey_list.append(np.nan)
            formula_list.append(np.nan)
            mw_list.append(np.nan)
            continue
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Row {idx+2}: Failed to parse SMILES: {smiles[:100]}...")
                clean_smiles_list.append(np.nan)
                parent_smiles_list.append(np.nan)
                inchi_list.append(np.nan)
                inchikey_list.append(np.nan)
                formula_list.append(np.nan)
                mw_list.append(np.nan)
                continue
            
            # 1. Remove salts (keep parent fragment)
            frag_parent = rdMolStandardize.FragmentParent(mol)
            clean_smiles = Chem.MolToSmiles(frag_parent, canonical=True, isomericSmiles=True)
            clean_smiles_list.append(clean_smiles)
            
            # 2. Neutralize charges to get parent structure
            uncharger = rdMolStandardize.Uncharger()
            uncharged_mol = uncharger.uncharge(frag_parent)
            parent_smiles = Chem.MolToSmiles(uncharged_mol, canonical=True, isomericSmiles=True)
            parent_smiles_list.append(parent_smiles)
            
            # 3. Generate InChI from salt-free molecule (some structures may fail)
            try:
                inchi = Chem.MolToInchi(frag_parent)
                # 4. Generate InChIKey
                if inchi:
                    inchikey = Chem.InchiToInchiKey(inchi)
                else:
                    inchikey = np.nan
            except Exception as e:
                print(f"Row {idx+2}: InChI generation failed: {e}")
                inchi = np.nan
                inchikey = np.nan
            inchi_list.append(inchi)
            inchikey_list.append(inchikey)
            
            # 5. Molecular formula
            formula = Chem.rdMolDescriptors.CalcMolFormula(frag_parent)
            formula_list.append(formula)
            
            # 6. Exact molecular weight
            mw = Descriptors.ExactMolWt(frag_parent)
            mw_list.append(round(mw, 4))
            
            # Progress indicator
            if (idx + 1) % 50 == 0:
                print(f"Processed {idx+1} compounds...")
                
        except Exception as e:
            print(f"Row {idx+2} error: {e}")
            clean_smiles_list.append(np.nan)
            parent_smiles_list.append(np.nan)
            inchi_list.append(np.nan)
            inchikey_list.append(np.nan)
            formula_list.append(np.nan)
            mw_list.append(np.nan)
    
    # Update DataFrame
    df['SMILES_RDKit'] = clean_smiles_list
    df['SMILES_Parent_RDKit'] = parent_smiles_list
    df['InChI_RDKit'] = inchi_list
    df['InChIKey_RDKit'] = inchikey_list
    df['Molecular_Formula_RDKit'] = formula_list
    df['Molecular_Weight_RDKit'] = mw_list
    
    # Save result
    df.to_excel(output_file, index=False)
    print(f"\nProcessing complete. Results saved to: {output_file}")
    
    # Statistics
    total = len(df)
    success_count = sum(1 for x in clean_smiles_list if not pd.isna(x))
    print(f"Total compounds: {total}")
    print(f"Successfully processed: {success_count}")
    print(f"Failed: {total - success_count}")
    print(f"Rows with both SMILES columns empty: {both_empty_count}")
    
    return df

if __name__ == "__main__":
    # Resolve paths relative to this script: RDkit/scripts -> RDkit/input, RDkit/output
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

    input_file = os.path.join(PROJECT_ROOT, "input", "data_table_Biologic_core.xlsx")
    output_file = os.path.join(PROJECT_ROOT, "output", "data_table_Biologic_core_standardized.xlsx")

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print("=" * 60)
    print("Starting SMILES processing (auto-select SMILES_PubChem or SMILES_Original)...")
    print("=" * 60)
    
    df_result = process_smiles_data(input_file, output_file)
