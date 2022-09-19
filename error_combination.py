# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:47:33 2022

@author: ZR48SA
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()


import random
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from pathlib import Path

import Bio
from Bio import SeqIO


#%%
# De novo sequencing error type	        Novor           PEAKS         PepNovo
#                                  CID (%)	HCD (%)	CID (%)	HCD (%)	CID (%)	HCD (%)
# Substitution of 1 by 1 or 2 AAs 	4.9 	7.7 	2.9 	6.9 	8.5 	7.1 
# Substitution of 2 by 2 AAs 	   10.5 	10.9 	9.4 	12.1 	22.7 	16.8 
# Substitution of 3 by 3 AAs 	    5.4 	5.5 	4.7 	6.0 	6.4 	9.6 
# Inversion of 2 or 3 AAs 	       11.1 	10.0 	9.1 	13.5 	31.5 	21.6 
# Substitution of 2 by 3 AAs 	    3.2 	3.9 	2.5 	4.2 	4.2 	4.4 
# Substitution of 4 by 4 AAs 	    8.2 	11.8 	7.2 	9.4 	9.3 	12.2 
# Substitution of 5 by 5 AAs 	    7.2 	9.5 	7.3 	8.1 	4.7 	8.7 
# Substitution of 6 by 6 AAs 	    6.9 	8.6 	6.5 	7.4 	4.2 	6.9 
# Other (AA accuracy >50%) 	        3.3 	2.7 	1.8 	2.8 	1.2 	0.4 
# Other (AA accuracy 25–50%) 	    31.0 	23.9 	27.9 	22.5 	5.2 	8.3 
# Other (AA accuracy <25%) 	         8.3 	5.3 	20.7 	7.0 	2.0 	4.0 
#Briefings in Bioinformatics, Volume 19, Issue 5, September 2018, Pages 954–970, https://doi.org/10.1093/bib/bbx033

#create custom combined error sets according to median


pdf=pd.read_excel("error_chances.xlsx",engine='openpyxl')
mean_counts=pdf.iloc[3:,1:].astype(float).mean(axis=1).round(1)

files=[

"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s6_6c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s5_5c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s4_4c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_3c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_3c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_1c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s1_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c5r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c5r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c5r1.fasta"]

#%% 
for rep in range(1,4):
    combined=[]
    # Substitution of 1 by 1 or 2 AAs
    ls=[]
    for file in ["C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_1c25r1.fasta",
              "C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s1_2c25r1.fasta"]:
    
        with open(file,"r") as f:
            records=SeqIO.parse(f,"fasta")
            df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
            r=df.sample(63)
            r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
            ls.append(r)
    combined.append(pd.concat(ls).sample(63))
    
    #Substitution of 2 by 2 AAs
    file="C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_2c25r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(137)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r)
        
    #Substitution of 3 by 3 AAs
    file="C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_3c25r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(63)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r)
        
    # Inversion of 2 or 3 AAs 
    ls=[]
    for file in ["C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c5r1.fasta",
               "C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c5r1.fasta"]:
    
        with open(file,"r") as f:
            records=SeqIO.parse(f,"fasta")
            df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
            r=df.sample(161)
            r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
            ls.append(r)
    combined.append(pd.concat(ls).sample(161))
            
    # Substitution of 2 by 3 AAs 
    ls=[]
    for file in ["C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_2c25r1.fasta",
               "C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_3c25r1.fasta"]:
    
        with open(file,"r") as f:
            records=SeqIO.parse(f,"fasta")
            df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
            r=df.sample(37)
            r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
            ls.append(r)
    combined.append(pd.concat(ls).sample(37))
             
    # Substitution of 4 by 4 AAs
    file="C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s4_4c25r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(97)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r)
        
    # Substitution of 5 by 5 AAs
    file="C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s5_5c25r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(76)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r)

    # Substitution of 6 by 6 AAs    
    file="C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s6_6c25r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(68)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r)
        
    # Other (AA accuracy >50%), Other (AA accuracy 25–50%), Other (AA accuracy <25%) (use mutations for this) 	
    file="C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c5r1.fasta"
    with open(file,"r") as f:
        records=SeqIO.parse(file,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
        r=df.sample(298)
        r.index=[Path(file).stem.strip("1").strip("r")]*len(r)
        combined.append(r) #added one extra to add up to total of 1000
        
    c=pd.concat(combined)
    c.to_csv("combined_"+"r"+str(rep)+".tsv",sep="\t")
    
    with open ("combined_"+"r"+str(rep)+".fasta","w")as f:
        f.write("\n".join(">"+c["original"]+"\n"+c["simulated"])+"\n")
