#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 21:00:10 2026

@author: courtneybrown
"""

import pandas as pd 
import os
import streamlit as st

df = pd.read_csv('molecule_analysis.csv')

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


lab_pic = 'IMG_5254 (1).png'
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
    
image_folder = 'mol_images'

for i, (index, row) in enumerate(filtered_df.iterrows()):

        with cols[i % 3]:
        img_path = os.path.join(image_folder, f"{row['Name']}.png")
        
        st.markdown(f"<h3 style='color: purple; '>{row['Name']}</h3>", unsafe_allow_html=True)
        st.write(f"**This molecule is involved in:** {row['Category']}")
        
      
        st.image(img_path, use_container_width=True)
        
        st.write(f"**Formula:** {row['Formula']}")
        st.write(f"**Molecular Weight:** {row['MW']}")
        st.write(f"**Hydrogen Bond Donors:** {row['HBD']}")
        st.write(f"**Hydrogen Bond Acceptors:** {row['HBA']}")
















