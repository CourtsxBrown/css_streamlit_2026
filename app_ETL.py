#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 20:49:58 2026

@author: courtneybrown
"""

import pandas as pd 
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# EXTRACT 

def extract_mols(folder_path):
    mol_list = []
    filenames = sorted([f for f in os.listdir(folder_path) if f.endswith(".mol")])

    for filename in filenames:
        file_path = os.path.join(folder_path, filename)
        mol = Chem.MolFromMolFile(file_path)
        
        if mol:
            base_name = os.path.splitext(filename)[0].capitalize()
            mol.SetProp("_Name", base_name) 
            mol_list.append(mol)
        
    print(f"Extraction complete. Loaded {len(mol_list)} molecules.")
    return mol_list

mol_list = extract_mols(r'/Users/courtneybrown/Library/Mobile Documents/com~apple~CloudDocs/PhD/CHPC CSS/Streamlit/my_app')

# TRANSFORM

category_map = {
    "Adrenaline": "Productivity",
    "Caffeine": "Energy",
    "Cortisol": "Productivity",
    "Glucose": "Energy",
    "Dopamine": "Happiness",
    "Serotonin": "Happiness"
}

def transform_mols(mol_list):
    
    mol_data = []
    
    for mol in mol_list:
        
        try:
            name = mol.GetProp("_Name") 
        except:
            name = "Unknown"
    
        category = category_map.get(name, "Uncategorized")
        mol.SetProp("Category", category)
    
        mol_data.append({
            "Name": name,
            "Category": category,
            "Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "MW": round(Descriptors.MolWt(mol), 2),
            "HBD": Chem.rdMolDescriptors.CalcNumHBD(mol),
            "HBA": Chem.rdMolDescriptors.CalcNumHBA(mol),
            })
    
    print(f"Transformation complete. Made descriptives for {len(mol_data)} molecules.")
    
    return mol_data 

mol_data = transform_mols(mol_list)

# LOAD

def load_mols(mol_data):

    df = pd.DataFrame(mol_data)
    
    print("Final DataFrame loaded and ready for Streamlit.")
    return df

df = load_mols(mol_data)

df.to_csv("molecule_analysis.csv", index=False)

def save_mol_images(mol_list, output_folder="mol_images"):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created folder: {output_folder}")

    for mol in mol_list:
        if mol:

            name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
            
   
            file_path = os.path.join(output_folder, f"{name}.png")
            

            Draw.MolToFile(mol, file_path, size=(300, 300))
            
    print(f"Images saved successfully to '{output_folder}'.")

images = save_mol_images(mol_list, output_folder="mol_images")
















