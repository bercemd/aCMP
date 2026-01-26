#Note: This script requires MDanalysis python module
# B Dutagaci, December, 2025 - Modified from the Original script from MDAnalysis


import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.analysis.distances import distance_array

import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.dcd',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--selseg',type=str, default='P011=P012',
    help='First selection.')
parser.add_argument(
    '--out_path',type=str, default='rmsf.dat',
    help='Output path.')
arg = parser.parse_args()

def contacts_mda(u,selection):

    sellist = []
    for i in range(len(selection.split("="))):
        sellist.append(selection.split("=")[i])

    sel = u.select_atoms("segid %s"%sellist[0], updating=True)
    ref = u.select_atoms("segid %s"%sellist[1], updating=True)
    sel_res = sel.residues
    ref_res = ref.residues
    #empty array to stack all frame's counts of contacts at positions i,j
    count_stack = np.zeros((len(sel_res), len(ref_res)))

    for frame in u.trajectory:
        print(frame)
        
        k=1
        #iterating through each of the selection pairs
        for i in range(len(sellist)):
            for j in range(len(sellist)):
                if j>i:
                    #initial selection of desired atoms both whole subunit and atoms close to target comparison
                    sel = u.select_atoms("segid %s"%sellist[i], updating=True)
                    ref = u.select_atoms("segid %s"%sellist[j], updating=True)
                    
                    #getting residue names and residue number from selected groups
                    sel_res = sel.residues
                    sel_seg = sel_res.segids
                    sel_resids = sel_res.resids
                    sel_nums = sel_res.resnums
                    
                    ref_res = ref.residues
                    ref_seg = ref_res.segids
                    ref_resids = ref_res.resids
                    ref_nums = ref_res.resnums
                    
        
                    #array of 100's to be replaced at positions i,j by the minimum distance between residue i and residue j
                    res_dist = np.full((len(sel_res), len(ref_res)), 100)
        
                    #iterating through each of the contacting residues
                    for n in sel_res:
                      for m in ref_res:
                
                        #selecting the atoms of sel and ref residues
                        sel_atoms = n.atoms.positions
                        ref_atoms = m.atoms.positions
                
                        #calculating the minimum distance between the residues
                        distance = distance_array(sel_atoms, ref_atoms, backend='OpenMP')
                        min_dist = np.min(distance)
                        #replacing the minimum distance value of residue i and residue j at position [i,j] in res_dist
                        sel_i = (np.where(n == sel_res)[0][0])
                        ref_j = (np.where(m == ref_res)[0][0])
                        res_dist[sel_i, ref_j] = min_dist
         
                    #Calculating where in res_dist is less than or equal to a cutoff distance
                    cutoff=5
                    contacts = contact_matrix(res_dist, cutoff)
                    #Converting True to 1, and False to 0
                    contact_int = contacts * 1
                    #only includes if the pair has any contact
                    if np.any(contact_int != 0):
                      #adding the contacts for the pair in the frame to the total contacts in count_stack
                      count_stack = count_stack + contact_int
                      k=k+1
                
    #Averaging the total number of contacts by the number of frames in the trajectory and over the number of pairs that have any contact
    count_avg = count_stack / len(u.trajectory) / k
    #counts = count_avg.transpose()
    counts = count_avg
    return counts, sel_resids, ref_resids

dcds = arg.traj_path

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
contacts_arr, sela_resids, selb_resids = contacts_mda(U,arg.selseg)
contacts_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(contacts_arr)):
  for j in range(len(contacts_arr[i])):
    print(sela_resids[i],selb_resids[j],contacts_arr[i][j],file=outputfile)
outputfile.close()

