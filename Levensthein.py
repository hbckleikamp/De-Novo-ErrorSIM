# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:11:59 2022

@author: ZR48SA
"""



import random
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from pathlib import Path

#%% matplotlib color cycle

import seaborn as sns
palette=sns.color_palette("Paired",n_colors=10)
#palette=[tuple(list(p)+[0.7]) for p in palette] #add alpha

#%%

from Levenshtein import distance as levenshtein_distance

import Bio
from Bio import SeqIO

#%%

files=["C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s6_6c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s6_6c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s6_6c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s5_5c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s5_5c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s5_5c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s4_4c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s4_4c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s4_4c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_3c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_3c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_3c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_2c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_2c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s3_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_3c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_3c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_3c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_2c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_2c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_1c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_1c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s2_1c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s1_2c100r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s1_2c50r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/subsitute_s1_2c25r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c20r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c10r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c5r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w3c1r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c20r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c10r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c5r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/invert_w2c1r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c20r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c10r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c5r1.fasta",
"C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c1r1.fasta"]

medians=[] 


#get baseline distribution
file="C:/paper 4/new/SwissprotIL de novo error simulation/mutate_c5r1.fasta"
t=Path(file).stem.strip("1").strip("r")
with open(file,"r") as f:
    records=SeqIO.parse(f,"fasta")
    df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
df.original=df.original.apply(lambda x: x.split(";")[0])
df["levensthein"]=[levenshtein_distance(x[0],x[1]) for x in df.values]
df["length"]=df.original.apply(len)
df["l/l"]=df["levensthein"]/df["length"]
bdf=df #baseline
    





for file in files:
 


    with open(file,"r") as f:
        records=SeqIO.parse(f,"fasta")
        df=pd.DataFrame([[r.description,str(r.seq)] for r in records],columns=["original","simulated"])
    df.original=df.original.apply(lambda x: x.split(";")[0])
    df["levensthein"]=[levenshtein_distance(x[0],x[1]) for x in df.values]
    df["length"]=df.original.apply(len)
    df["l/l"]=df["levensthein"]/df["length"]

    
    fig,ax=plt.subplots()
    plt.hist(bdf["l/l"],color=palette[0]) #baseline (5 chance mutation)
    plt.hist(df["l/l"],color=list(palette[1])+[0.75]) # add alpha
    plt.xlabel("Levensthein distance / peptide length")
    plt.ylabel("Frequency")

    plt.title(Path(file).stem)
    plt.text(0.95, 0.5, 'mutate_c5 median='+str(round(bdf["l/l"].median(),3)), horizontalalignment='right',verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95, 0.4, Path(file).stem.strip("1").strip("r")+' median='+str(round(df["l/l"].median(),3)), horizontalalignment='right',verticalalignment='center', transform=ax.transAxes)
    plt.legend(labels=[t,Path(file).stem.strip("1").strip("r")])
    plt.savefig("Levensthein_"+Path(file).stem.strip("1").strip("r")+".png",dpi=600,bbox_inches="tight")

    medians.append([Path(file).stem,df["l/l"].median()])
    
mdf=pd.DataFrame(medians,columns=["file","median levensthein distance / peptide length"])
mdf["chance"]=mdf.file.apply(lambda x: x.split("c")[1].split("r")[0])
mdf.to_csv("median levensthein distance.csv")