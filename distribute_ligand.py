#! /usr/bin/env python3

import random
import sys
import pymol
import os


if len(sys.argv) < 7 :
    print("Usage : ./distribute_ligand pdb ligand_3_code nbr_copies angstrom_range target_name rm_flag")
    exit(1)

#give .pdb filename in first argument
prot_path = os.path.abspath(sys.argv[1])
prot_name = sys.argv[1].split('.')[0]
target_name = sys.argv[5]
#give ligand three letter code in second parameter
lig_code = sys.argv[2]
nbr_copies = sys.argv[3]
ang_rang = int(sys.argv[4])
#boolean flag to remove initial ligand or not 0 no 1 yes 
rm_flag = int(sys.argv[6])

#load pdb
pymol.cmd.load(prot_path, prot_name)

#cleanup water because new molecules will be added
pymol.cmd.remove("sol.")
pymol.cmd.remove("not chain A")

#center the protein
pymol.cmd.center("all")

#get box coordinates around protein
#using pymol
([min_X, min_Y, min_Z],[max_X, max_Y, max_Z]) = pymol.cmd.get_extent("all")

#copy, paste, translate and rotate ligand
pymol.cmd.select("lig_1", f"resname {lig_code}")
pymol.cmd.extract("ligand", "lig_1")
ligands = []
#need to generate range between protein box and !angstrom_range
i=0
if int(nbr_copies) > 0 :
    while i < int(nbr_copies) :
        j=i+2
        pymol.cmd.copy(f"lig{j}", "ligand")
        axe1 = random.randrange(0,3,1)
        axe2 = random.randrange(0,3,1)
        x_sign = random.getrandbits(1)
        y_sign = random.getrandbits(1)
        z_sign = random.getrandbits(1)
        
        #is X
        if axe1 == 0 :
            if x_sign :
                range_distrib_x = random.uniform(min_X - ang_rang, min_X-1)
                if axe2 == 1 :
                    if y_sign :
                        range_distrib_y = random.uniform(min_Y - ang_rang, min_Y)
                    else :
                        range_distrib_y = random.uniform(max_Y, max_Y + ang_rang)
                else :
                    range_distrib_y = random.uniform(min_Y, max_Y)
            else :
                range_distrib_x = random.uniform(max_X+1, max_X + ang_rang)
                if axe2 == 1 :
                    if y_sign :
                        range_distrib_y = random.uniform(min_Y - ang_rang, min_Y)
                    else :
                        range_distrib_y = random.uniform(max_Y, max_Y + ang_rang)
                else :
                    range_distrib_y = random.uniform(min_Y, max_Y)
            range_distrib_z = random.uniform(min_Z, max_Z)

        #is Y
        if axe1 == 1 :
            if y_sign :
                range_distrib_y = random.uniform(min_Y - ang_rang, min_Y-1)
                if axe2 == 1 :
                    if z_sign :
                        range_distrib_z = random.uniform(min_Z - ang_rang, min_Z)
                    else :
                        range_distrib_z = random.uniform(max_Z, max_Z + ang_rang)
                else :
                    range_distrib_z = random.uniform(min_Z, max_Z)
            else :
                range_distrib_y = random.uniform(max_Y+1, max_Y + ang_rang)
                if axe2 == 1 :
                    if z_sign :
                        range_distrib_z = random.uniform(min_Z - ang_rang, min_Z)
                    else :
                        range_distrib_z = random.uniform(max_Z, max_Z + ang_rang)
                else :
                    range_distrib_z = random.uniform(min_Z, max_Z)
            range_distrib_x = random.uniform(min_X, max_X)

        
        #is Z
        if axe1 == 2 :
            if z_sign :
                range_distrib_z = random.uniform(min_Z - ang_rang, min_Z-1)
                if axe2 == 1 :
                    if x_sign :
                        range_distrib_x = random.uniform(min_X - ang_rang, min_X)
                    else :
                        range_distrib_x = random.uniform(max_X, max_X + ang_rang)
                else :
                    range_distrib_x = random.uniform(min_X, max_X)
            else :
                range_distrib_z = random.uniform(max_Z+1, max_Z + ang_rang)
                if axe2 == 1 :
                    if x_sign :
                        range_distrib_x = random.uniform(min_X - ang_rang, min_X)
                    else :
                        range_distrib_x = random.uniform(max_X, max_X + ang_rang)
                else :
                    range_distrib_x = random.uniform(min_X, max_X)
            range_distrib_y = random.uniform(min_Y, max_Y)

        
        #would need to note all the positions to verify not two ligands are at the same place
        #maybe while instead of for with a control of distances?
        coords = [range_distrib_x,range_distrib_y,range_distrib_z] 
        ligands.append(coords)
        
        
        temp = 0
        for k in range(len(ligands)-1) :
                diff = abs(ligands[k][0]-coords[0])+abs(ligands[k][1]-coords[1])+abs(ligands[k][2]-coords[2])
                temp = max(temp,diff)

        if temp < 4 : # constant? or to determine
            #skip a loop
            continue

        #only if good
        i+=1
        pymol.cmd.translate(coords,f"lig{j}")
        

        pymol.cmd.rotate("z",random.randrange(0,360,1),f"lig{j}")
        pymol.cmd.rotate("y",random.randrange(0,360,1),f"lig{j}")
        pymol.cmd.rotate("x",random.randrange(0,360,1),f"lig{j}")


#this is good enough, but all ligands will be considered
#a single entity except if the file is manually altered
if rm_flag :
    pymol.cmd.remove("ligand")
pymol.cmd.save("temp.pdb", "all", -1, "pdb")

#here is a try for the cleanup
#count number of atoms in ligand
switch=0
i=0
with open(sys.argv[1],"r") as f:
    for line in f :
        line = line.rstrip("\n")
        line = line.split()
        if len(line) > 3 :
            if line[3] == lig_code :
                switch=1
                i+=1
            if switch == 1:
                if line[3] != lig_code :
                    break
stock = i
ligg = int(nbr_copies) if 1 else int(nbr_copies + 1)
#to separate ligands
with open(f"{target_name}.pdb", "w+") as fo :
    with open("temp.pdb", "r") as fi :
        for line in fi:
            linee = line.rstrip("\n")
            linee = linee.split()
            if len(linee) > 3 :
                if linee[3] == lig_code:
                    j=9999
                    if linee[0] != "CONECT" :
                        
                        if linee[3] == lig_code :
                            #sorry
                            blipou = str(j-ligg)
                            line = list(line)
                            line[22:26] = [k for k in blipou]
                            line = ''.join(line)
                            
                            i-=1
                            fo.write(line)
                            if i ==  0 :               
                                ligg -= 1
                                i = stock
                                fo.write("TER\n")
                    
                else :
                    fo.write(line)
                    

            else :
                fo.write(line)


pymol.cmd.quit()