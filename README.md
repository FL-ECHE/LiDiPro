# LiDiPro
Ligand Distributor around Protein

Script using Pymol to distribute ligands of interest around a protein randomly.

## Requirements

Pymol must be installed on the computer, I use it as a python library to first select and then translate and rotate copies of the ligand.

You should know the three letter code figuring in the PDB for the ligand of choice, for it is one of the inputs

## Use

Call the script on the pdb file with ligand PDB code, number of copies desired, range around the protein where the ligands will be distributed, name of the target file where the new complex will be saved, and wether to delete the original ligand (0 for no, 1 for yes).

This intended to easily set up a system for binding molecular dynamics simulations.
