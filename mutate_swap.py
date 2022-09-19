# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 13:54:55 2022

@author: ZR48SA
"""



import random
import numpy as np
import pandas as pd
import re



# Evaluating de novo sequencing in proteomics: already an accurate alternative to database-driven peptide identification? 
# Thilo Muth, Bernhard Y Renard
# Briefings in Bioinformatics, Volume 19, Issue 5, September 2018, Pages 954–970, https://doi.org/10.1093/bib/bbx033


#type                              #average chance
# Substitution of 1 by 1 or 2 AAs 6.3
# Substitution of 2 by 2 AAs    13.7
# Substitution of 3 by 3 AAs      6.3
# Inversion of 2 or 3 AAs        16.1
# Substitution of 2 by 3 AAs        3.7
# Substitution of 4 by 4 AAs     9.7
# Substitution of 5 by 5 AAs      7.6
# Substitution of 6 by 6 AAs       6.8
# Other (AA accuracy >50%)         2.0
# Other (AA accuracy 25–50%)    19.8
# Other (AA accuracy <25%)         7.9


#%% general
std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'N': 114.04293,
               'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,'M': 131.04049,
               'H': 137.05891,'F': 147.06841,'R': 156.10111,'Y': 163.06333,
               'W': 186.07931}    

def pep_mass_calc(x,std_aa_mass=std_aa_mass):
    return sum(std_aa_mass.get(aa) for aa in x)+18.01056

c2a = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

a2c = {y:x for x,y in c2a.items()}

#%% Mutation

def mutate(peptide,chance=5,c2a=c2a,a2c=a2c):

    #reverse translation
    rt="".join([a2c.get(p) for p in peptide]) 
    
    #mutation
    chance=int(chance)
    mrt="".join(["ACGT"[random.randint(0,3)] if random.randrange(1, 101, 1)<=chance*4/3 else rt[ix] for ix,i in enumerate(range(len(rt)))])
    
    #translation
    tr="".join([c2a.get(mrt[i:i+3]) for i in range(0, len(mrt), 3)]) 
    return tr


#%% Mass substitution

import itertools
a = [i for i in std_aa_mass.keys()]
r=[]
for i in range(1,7):
    r.extend(list(itertools.combinations_with_replacement(a,i)))
masses=[sum([std_aa_mass.get(aa) for aa in i]) for i in r]
peptides=["".join(sorted(i)) for i in r]        
subdf=pd.DataFrame(list(zip(masses,peptides)),columns=["mass","peptide"])
subdf["mass"]=subdf["mass"].round(3)
subdf["peptide"]=subdf["peptide"].apply(lambda x: x.replace("I","L"))
subdf=subdf.drop_duplicates()
subs=subdf.groupby("mass").apply(len)
subdf=subdf[subdf["mass"].isin(subs[subs>1].index)].sort_values(by="mass")
subdf=subdf.set_index("mass")
subdf["l"]=subdf["peptide"].apply(len)



def substitute(peptide,subno,nosub,chance,subdf=subdf):
     
    #find subsitutable locations
    slide=["".join(peptide[i:i+subno]) for i in range(len(peptide)-(subno-1))]
    sd=subdf.loc[subdf[subdf.peptide.isin(slide)].index,:]
    targets=sd[(~sd.peptide.isin(slide)) & (sd.l==nosub)]
    s1=sd.loc[targets.index,:]
    s2=s1[s1.peptide.isin(slide)]
    s3=s2.peptide.tolist()
    sx=[ix for ix,i in enumerate(slide) if i in s3]
    
    
    #pick substitutable locations that do not overlap, with supplied chance 
    ixs=[]
    iy=0
    for ix in sx:
        if ix>=iy:
            if random.randint(0, 100)<chance:
                ixs.append(ix)
                iy=ix+subno
        
        if ix>len(slide):
            break
    
    #substitute peptide at non-overlapping locations with randomized equal mass substitution
    iy=0
    pm=""
    for ix,i in enumerate(peptide):
        if ix>=iy:
            if ix in ixs:
                s=list(targets.loc[s2[s2.peptide==slide[ix]].index].sample(1).peptide.tolist()[0])
                random.shuffle(s)
                s="".join(s)
                pm+=s
                iy=ix+subno
            else:    
                pm+=i[0]
    return pm

#%% inversion

def invert(peptide,subno,chance):

    slide=["".join(peptide[i:i+subno]) for i in range(len(peptide)-(subno-1))]
    #pick substitutable locations that do not overlap, with supplied chance 
    ixs=[]
    iy=0
    for ix,i in enumerate(slide):
        if ix>=iy:
            if random.randint(0, 100)<chance:
                ixs.append(ix)
                iy=ix+subno
        if ix>len(slide):
            break
    #substitute peptide at non-overlapping locations with randomized equal mass substitution
    iy=0
    pm=""
    for ix,i in enumerate(peptide):
        if ix>=iy:
            if ix in ixs:
                s=list(slide[ix])
                random.shuffle(s)
                pm+="".join(s)
                iy=ix+subno
            else:    
                pm+=i[0]
    return pm



#%% 
#loop over peptides in dataset


file="C:/paper 4/new/SwissprotIL precomutping LCA/ncbi/swisprotIL_precomputed_2.tsv"
pdf=pd.read_csv(file,sep="\t")

for a in ["B","Z","J","X","O","U"]: #remove ambiguous amino acids and unusual amino acids B,Z,J,X,U,O
    pdf=pdf[~pdf.peptides.str.contains(a)]


#%%
for rep in range(1,2): #range(1,4)
    
    subset=pdf.sample(100000).fillna("")
    peptides=subset.peptides
    header=subset.apply(";".join,axis=1)

    #mutate
    for chance in [1,5,10,20]:
        p=peptides.apply(mutate,chance=chance)
        o=p[~p.isin(peptides)].sample(1000)
        
        with open ("mutate_c"+str(chance)+"r"+str(rep)+".fasta","w")as f:
            f.write("\n".join(">"+header[o.index]+"\n"+o)+"\n")

    #invert
    for window in [2,3]:
        for chance in[ 1,5,10,20]:
            p=peptides.apply(invert,subno=window,chance=chance)
            o=p[~p.isin(peptides)].sample(1000)
            
            with open ("invert_w"+str(window)+"c"+str(chance)+"r"+str(rep)+".fasta","w")as f:
                f.write("\n".join(">"+header[o.index]+"\n"+o)+"\n")

    #invert
    for subs in [[1,2],
                 [2,1],
                 [2,2],
                 [2,3],
                 [3,2],
                 [3,3],
                 [4,4],
                 [5,5],
                 [6,6]]:
        
        
        sdf=subdf[(subdf.l==subs[0]) | (subdf.l==subs[1])]
        
        for chance in[25,50,100]:
            p=peptides.apply(substitute,subno=subs[0],nosub=subs[1],chance=chance,subdf=sdf)
                        
            o=p[~p.isin(peptides)].sample(1000)
            
            with open ("subsitute_s"+str(subs[0])+"_"+str(subs[1])+"c"+str(chance)+"r"+str(rep)+".fasta","w")as f:
                f.write("\n".join(">"+header[o.index]+"\n"+o)+"\n")

  
        

