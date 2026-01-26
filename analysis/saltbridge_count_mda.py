#Note: This script requires MDanalysis python module
# B Dutagaci, December, 2025 


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
    '--sela',type=str, default='resname ARG and (name NH1 or name NH2)',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='resname GLU and (name OE1 or name OE2)',
    help='Second selection.')
parser.add_argument(
    '--out_path',type=str, default='salt_bridges.dat',
    help='Output path.')
arg = parser.parse_args()

def counts_mda(u,selecta,selectb):
    sela = u.select_atoms(selecta)
    selb = u.select_atoms(selectb)
    n_sela = len(sela)
    n_selb = len(selb)

    print('Sela has {} residues and Selb has {} residues'.format(n_sela, n_selb))

    sela_res = sela.residues
    selb_res = selb.residues
    sela_resids = sela_res.resids
    selb_resids = selb_res.resids
    count_arr = np.zeros(shape=(len(u.trajectory)))
    boxsize = u.trajectory.ts.dimensions
    i=0
    for frame in u.trajectory:
      print(frame)
      count = 0
      distance = distance_array(sela.atoms.positions,selb.atoms.positions,box=boxsize,backend='OpenMP')
      cutoff=4
      contacts = contact_matrix(distance, cutoff)
      contact_int = contacts * 1
      count += int(np.sum(contact_int))

      count_arr[i] = count
      i=i+1

    return count_arr

dcds = arg.traj_path

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
count_arr = counts_mda(U,arg.sela,arg.selb)
count_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(count_arr)):
    print(count_arr[i],file=outputfile)
outputfile.close()

