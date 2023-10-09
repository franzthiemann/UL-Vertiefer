#!/usr/bin/python

# Uses geodesic interpolation for the molecule and idpp interpolation for the surface of a molecule

from ase.io import read, write
from ase import Atoms
import re
import os
import sys
from ase.neb import NEB
from ase.calculators.lj import LennardJones as LJ
from ase.io import read,Trajectory
from ase.visualize import view
from pathlib import Path
import csv
from decimal import Decimal
import ase.io
from geodesic_interpolate.interpolation import redistribute
from geodesic_interpolate.geodesic import Geodesic

######## Functions as copied from geodesic_wrapper.py ########
def ase_geodesic_interpolate(initial_mol,final_mol, n_images = 20, friction = 0.01, dist_cutoff = 3, scaling = 1.7, sweep = None, tol = 0.002, maxiter = 15, microiter = 20):
    atom_string = initial_mol.symbols
    
    atoms = list(atom_string)

    initial_pos = [initial_mol.positions]
    final_pos = [final_mol.positions]

    total_pos = initial_pos + final_pos

    # First redistribute number of images.  Perform interpolation if too few and subsampling if too many
    # images are given
    raw = redistribute(atoms, total_pos, n_images, tol=tol * 5)

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    smoother = Geodesic(atoms, raw, scaling, threshold=dist_cutoff, friction=friction)

    if sweep is None:
        sweep = len(atoms) > 35
    try:
        if sweep:
            smoother.sweep(tol=tol, max_iter=maxiter, micro_iter=microiter)
        else:
            smoother.smooth(tol=tol, max_iter=maxiter)
    finally:

        all_mols = []
        
        for pos in smoother.path:
            mol = Atoms(atom_string, pos)
            all_mols.append(mol)
        
        return all_mols

########## Functions copied from geodesic_wrapper.py ##########
def interpolate_traj(initial,final,LEN,calculator=None):
    interpolated_length=LEN
    initial_1=initial.copy()
    initial_1.calc=calculator #attach calculator to images between start and end
    RETURN=initial_1.copy()
    images = [initial]
    images += [initial_1.copy() for i in range(interpolated_length)] #create LEN instances
    images += [final]
    nebts = NEB(images)
    nebts.interpolate(mic=True,method='idpp',apply_constraint=True)
    #nebts.idpp_interpolate(mic=True)
    new_traj=Trajectory(f'{LEN}_interpol.traj',mode='w')
    for im in images:
        new_traj.write(im)
def create_trajs(START,END,LEN):
    calc = LJ()
    #TRAJ=interpolate_traj(START,END,30,calculator=calc)
    #TRAJ=interpolate_traj(START,END,20,calculator=calc)
#    TRAJ=interpolate_traj(START,END,12,calculator=calc)
    TRAJ=interpolate_traj(START,END,LEN,calculator=calc)

######### New Functions #########

def main():
    if len(sys.argv) < 5:
        print("Uses geodesic interpolation for the molecule and idpp interpolation for the surface of a molecule")
        print("Usage:")
        print("interpolateMixed.py [start] [end] [Images] [SurfaceCutoff]")
        print("[start]:             POSCAR file of initial molecule")
        print("[end]:               POSCAR file of final molecule")
        print("[Images]:            Number of NEB Images to generate")
        print("[SurfaceCutoff]:     Atom Number of the first Atom that belongs to the molecule")
        sys.exit()
    initial_mol=read(sys.argv[1])
    atoms=initial_mol
    final_mol=read(sys.argv[2])
    LEN=int(sys.argv[3])
    moleculeStart = int(sys.argv[4])
    if len(sys.argv) == 6 and sys.argv[5] == "True":
        print("using alternate mode")
        initial_surf = initial_mol[:moleculeStart]
        final_surf = initial_mol[:moleculeStart]
        create_trajs(initial_surf,final_surf,LEN)
    else:
        create_trajs(initial_mol,final_mol,LEN)
    trajectory_idpp=read(f'{LEN}_interpol.traj',index=':')
    trajectory_geodisic = ase_geodesic_interpolate(initial_mol,final_mol, n_images= LEN+1)
    difference=atoms[0].position-trajectory_geodisic[0][0].position
    for image in trajectory_geodisic:
        image.set_cell(atoms.cell)
        for at in image:
            at.position=at.position+difference
    # For some reason, the geosidic trajectory misses the last frame. Increasing the number of frames does not help since the last
    # frame is not the POSCAR.end
    newTrajectory = trajectory_idpp.copy()
    trajectory_geodisic.append(final_mol)
    for nr, image in enumerate(trajectory_idpp):
        atoms_idpp = trajectory_idpp[nr]
        atoms_geodisic = trajectory_geodisic[nr]
        molecule = atoms_geodisic[moleculeStart:]
        surface = atoms_idpp[:moleculeStart]
        newTrajectory[nr] = Atoms(surface + molecule)
        merged = Atoms(surface + molecule)
        N=f'{nr:02d}'
        Path(N).mkdir(parents=True,exist_ok=True)
        newTrajectory[nr].write(f'{N}/POSCAR',format='vasp')
main()
