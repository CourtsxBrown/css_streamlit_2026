#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 16:28:40 2026

@author: courtneybrown
"""

import pandas as pd 
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

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
            "Object": mol,
            "SMILES": Chem.MolToSmiles(mol)
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

import streamlit as st 
from rdkit.Chem import Draw

st.markdown(
 """
<style>
 /* Targets only the main titles */
 h1 {
     color: black;
 }
/* Targets only the subheaders */
h3 {
    color: purple!important;
}
 .stApp {
     background-color: white;
     color: black
     }
[data-testid="stSidebar"] {
            background-color: grey;
        }
 </style>
 """,
 unsafe_allow_html=True
 )


lab_pic = r'/Users/courtneybrown/Library/Mobile Documents/com~apple~CloudDocs/PhD/CHPC CSS/Streamlit/my_app/IMG_5254 (1).png'

col1, col2 = st.columns([3, 5])
with col1:
    st.image(lab_pic, width=300)
with col2: 
    st.title("BIOCHEMIST | MSc Candidate")
    st.write("Courtney Brown")
    st.markdown('<p style="font-size: 20px; "><b>Experience life as a Biochemist</b> — <i>Explore the chemical properties of a compound library.</i></p>', unsafe_allow_html=True)
    
st.sidebar.header("Molecule Filtering Criteria")
    
categories = list(df["Category"].unique())

available_mols = df["Name"].tolist()
    
selected_mol = st.sidebar.multiselect("Select Specific Molecules", available_mols)

mw_limit = st.sidebar.number_input("Maximum Molecular Weight", value=600.0, step=50.0)

filtered_df = df[df["Name"].isin(selected_mol)]
filtered_df = filtered_df[filtered_df["MW"] <= mw_limit]

if filtered_df.empty:
        st.markdown("**Choose your filtering criteria using the side pannel.**")
        st.markdown("(Tip: You can look at more than one molecule at a time.)")
else:
    cols = st.columns(3)

for index, row in filtered_df.iterrows():
    with cols[index % 3]:
        st.markdown(f"<h3 style='color: purple; '>{row['Name']}</h3>", unsafe_allow_html=True)
        st.write(f"**This molecule is involved in:** {row['Category']}")
        st.write(f"**Formula:** {row['Formula']}")
        st.write(f"**Molecular Weight:** {row['MW']}")
        st.write(f"**Hydrogen Bond Donors:** {row['HBD']}")
        st.write(f"**Hydrogen Bond Acceptors:** {row['HBA']}")
        
        smiles_str = row['SMILES']
        
        mol = Chem.MolFromSmiles(smiles_str)
        
        if mol:
            img = Draw.MolToImage(mol, size=(500, 500))
            st.image(img, use_container_width=True)
        else:
            st.warning(f"Could not render structure for {row['Name']}")
                
                
