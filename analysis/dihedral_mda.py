#!/usr/bin/env python3

#Note: This script requires MDanalysis python module
# B Dutagaci, 2025

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.dihedrals as dihedrals
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
    '--selseg',type=str, default='segid P011',
    help='First selection.')
parser.add_argument(
    '--out_path',type=str, default='rmsd.dat',
    help='Output path.')
arg = parser.parse_args()

u = mda.Universe(arg.psf_path, arg.traj_path)
selres = u.select_atoms(arg.selseg)

phi = [res.phi_selection() for res in selres.residues]
Rphi = dihedrals.Dihedral(phi).run()

psi = [res.psi_selection() for res in selres.residues]
Rpsi = dihedrals.Dihedral(psi).run()

outputfile = open(arg.out_path,"w")

for i in range(len(Rphi.angles)):
    for j in range(len(Rphi.angles[i])):
        print(Rphi.angles[i][j],Rpsi.angles[i][j],file=outputfile)

outputfile.close()

