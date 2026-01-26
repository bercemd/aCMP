#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# J Oakden, B Dutagaci, 2025 


import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")
from collections import Counter
from tqdm import tqdm
import pandas as pd
import os,sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
        '--title',type=str,default="title",
        help='title_error'
        )
parser.add_argument(
        '--psf',type=str,default="protein.psf",
        help='psf file'
        )
parser.add_argument(
        '--dcd',type=str,
        help='dcd file'
        )
parser.add_argument(
        '--selection',type=str,default="name CA",
        help='fit selection'
        )
parser.add_argument(
        '--distance',type=float,default=7.0,
        help='cut off'
        )
arg =parser.parse_args()

PSF =arg.psf
DCDpath =arg.dcd
Name =arg.title
cutoff=arg.distance
select=arg.selection
DCD=[]
lis=DCDpath.split(" ")
for i in lis:
    DCD.append(i)
u = mda.Universe(PSF, DCD)
boxsize = u.trajectory.ts.dimensions

prot = u.select_atoms(f"protein")
prot_seg = prot.segments
prot_id = prot_seg.segids
seglis=[]
for segment in prot_id:
    seglis.append(segment)
distime=[]
pairselection=[]
comb = list(combinations(seglis,2))
for i in comb:
    iseg=u.select_atoms(f"segid {i[0]} and {select}")
    jseg=u.select_atoms(f"segid {i[1]} and {select}")
    a=[i,iseg,jseg]
    pairselection.append(a)
DATALIST=[]
for ts in u.trajectory:
    for seg in pairselection:
        mydist=distances.distance_array(seg[1].atoms.positions, seg[2].atoms.positions, box = boxsize , backend="OpenMP")


        minval=np.min(mydist)
        if minval < cutoff:
            DATALIST.append(f"{ts.frame} {seg[0][0]} {seg[0][1]} {minval}")
        else:
            DATALIST.append(f"{ts.frame} {seg[0][0]} {seg[0][1]}")
numseg=len(prot_id)

amount=[0]*numseg
listofweird = []
normallist = []
last = None 
            
        
for line in tqdm(DATALIST):
    parts = line.strip().split()
    if not parts:
        continue
    last = int(parts[0])
    if len(parts) == 4:
        normallist.append(parts)
    elif len(parts) > 4:
        listofweird.append(parts)
if listofweird:
    print("uh oh")

frames = {}
for data in tqdm(normallist):
    index = data[0]
                
    if index not in frames:
        frames[index] = []
    frames[index].append(data)

clust = []
for index in tqdm(frames):
                
    framedat = frames[index]
    pairs = [[seg[1], seg[2]] for seg in framedat]
                
    clusters = []
    for a, b in pairs:

        match = []
                    
        indx = 0
        for cluster in clusters:
            if a in cluster or b in cluster:
                match.append(indx)
            indx=indx+1
                    
        if len(match) == 0:
            clusters.append({a,b})
        else:
            clean = {a,b}            
            match.sort()
            match.reverse()
            
            for tir in match:
                oldclus = clusters.pop(tir)
                clean.update(oldclus) 
            clusters.append(clean)

    ficlusters = [clus for clus in clusters]
    clust.append([index, ficlusters])

    count = 0
    for clus in ficlusters:
        nestlistlen = len(clus)
        amount[nestlistlen-1] = amount[nestlistlen-1]+1
        count=count+int(nestlistlen)      
    if count<numseg:
        amount[0]=amount[0]+(numseg-count)
gowithit = len(frames)
frwnoclus=(last + 1) - gowithit
amount[0]=amount[0] + (frwnoclus * numseg)
tot=0
for i in range(len(amount)):
  tot = tot+(amount[i]*(i+1))
percentage = []
for i in range(len(amount)):
  percentage.append((amount[i]*100*(i+1))/tot)

with open(f"{Name}.pcnt.txt", "w") as f:
    for i in range(len(amount)):
      j=i+1  
      f.write(f"{j} {amount[i]} {percentage[i]}\n")
 
