# LiDiPro
Ligand Distributor around Protein

Script using Pymol to distribute ligands of interest around a protein randomly.

## Requirements

Pymol must be installed on the computer, I use it as a python library to first select and then translate and rotate copies of the ligand.

You should know the three letter code figuring in the PDB for the ligand of choice, for it is one of the inputs

## Use

Call the script on the ***pdb file*** with ***ligand PDB code***, ***number of copies desired***, ***range around the protein*** where the ligands will be distributed, ***name of the target file*** where the new complex will be saved, and wether to delete the original ligand ***(0 for no, 1 for yes)***.

Example to distribute 5 copies of the XYZ ligand around 5 angström of the complex :
```bash
./distribute_ligand.py test.pdb XYZ 5 5 test_5lig.pdb 0 
```

This project is intended to easily set up a system for binding molecular dynamics simulations.

It allows to make copies of a specific molecule
![one ligand](lidipro1.png)

And randomly distribute them around the protein of interest
![several ligand](lidipro2.png)

For now for my own usage, I only need the A chain, so I remove the rest, but if you need the whole complexed subunits, you can suppress this line :

```python
pymol.cmd.remove("not chain A")
```
