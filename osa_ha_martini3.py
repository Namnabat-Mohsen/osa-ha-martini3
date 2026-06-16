# =============================================================================
# osa_ha_martini3.py
#
# Generator of primitive starting geometries and Martini 3 force-field files
# (GROMACS format) for hyaluronic acid (HA) and OSA-grafted HA (OSA-HA).
#
# -----------------------------------------------------------------------------
# ATTRIBUTION / LICENSE
# -----------------------------------------------------------------------------
# This file is a DERIVATIVE WORK of `carbo2martini3_2.0.py` by
#   Valery Lutsyk and Wojciech Plazinski,
#   Supporting Information of:
#   "Extending the Martini 3 Coarse-Grained Force Field to Hyaluronic Acid",
#   J. Phys. Chem. B 2025, 129 (9), 2408-2425.
#   DOI: 10.1021/acs.jpcb.4c08043
# The original work is distributed by the authors under CC BY 4.0
# (https://creativecommons.org/licenses/by/4.0/).
#
# The Martini 3 carbohydrate machinery this descends from was first described in:
#   V. Lutsyk, P. Wolski, W. Plazinski,
#   "Extending the Martini 3 Coarse-Grained Force Field to Carbohydrates",
#   J. Chem. Theory Comput. 2022, 18 (8), 5089-5107.
#   DOI: 10.1021/acs.jctc.2c00553
#
# CHANGES MADE relative to the original carbo2martini3_2.0.py:
#   - Added the OSA-grafted HA generator `OSA_HA(units_num, ds_percent, seed)`
#     and its supporting functions (`_ha_osa_templates`, `_ha_absolute_coord`,
#     `_place_osa_fragment`, `_resolve_repeat_local`).
#   - Added octenyl-succinate (OSA) bead definitions and the ester-linkage
#     bonded terms attaching the OSA fragment to the GlcNAc residue.
#   - Added degree-of-substitution (DS) control with a reproducible random seed,
#     and menu option 7 ("OSA_HA") in main().
#   The original HA / glucopyranose force-field parameters and the core
#   topology-writing machinery are UNCHANGED.
#
# Modifications (c) 2026 [YOUR FULL NAME], released under CC BY 4.0.
# This derivative is NOT endorsed by the original authors.
#
# >> SCIENTIFIC CAVEAT <<
# The OSA graft bead types and the ester-linkage bonded parameters are a
# PRELIMINARY, UNVALIDATED extension. Unlike the published HA parameters
# (validated against atomistic MD and octanol-water transfer free energies),
# the OSA additions here have NOT been independently validated. Validate
# against an atomistic reference before any production use.
# =============================================================================

import numpy as np
import math
import random

"""
A general definition of the single residue common to glucose-based polysaccharides:
     1
    /
   4 -- 2
    \  /
      3  

A general definition of the single residue for hyaluronic acid:
     1         6
    /         /
   4 -- 2 -- 9 -- 7
    \  /      \  / 
      3         8
       \
        5
"""



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def pickAngle(index):
    if (index+1) % 4 == 1:
        angle = 0
    elif (index+1) % 4 == 2:
        angle = -1.5707963268
    elif (index+1) % 4 == 3:
        angle = -3.1415926536
    elif (index+1) % 4 == 0:
        angle = -4.7123889804
    return angle


class Polysaccharide:
    """
    The class that stores the parameters of a molecule composed of multiple residues.
    The class also has write_gro and write_top methods to create gro and top files
    """

    def __init__(self,name, unit_at):

        self.name      = name
        
        self.atoms     = []
        self.coordsX   = []
        self.coordsY   = []
        self.coordsZ   = []
        self.bonds     = []
        self.angles    = []
        self.dihedrals = []
        self.impropers = []
        
        self.qtot      = 0
        self.unit_at   = unit_at

    def gets_atoms(self, index, row):
        self.qtot += row[4]
        self.atoms.append([row[0]+index*self.unit_at, row[1], index+1, row[2], row[3], row[0]+index*self.unit_at, row[4], row[5], self.qtot])

    def gets_bonds(self, index, row):
        self.bonds.append([row[0]+index*self.unit_at,row[1]+index*self.unit_at,row[2],row[3]])

    def gets_angles(self, index, row):
        self.angles.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3], row[4], row[5]])

    def gets_dihedrals(self, index, row):
        self.dihedrals.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3]+index*self.unit_at, row[4], row[5], row[6], row[7]])
    
    def gets_impropers(self, index, row):
        self.impropers.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3]+index*self.unit_at, row[4], row[5], row[6]])

    
    
    def writeGro(self):
        
        fname = self.name + ".gro"
        groFile = open(fname, "w")
        groFile.write("%s\n" % self.name)
        groFile.write(" %i\n" % len(self.atoms))
        count = 1

        for at, x, y, z  in zip(self.atoms,self.coordsX, self.coordsY, self.coordsZ):
            line = "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (at[2], at[3], at[4], count, x, y, z)
            count += 1
            if count == 100000:
                count = 0
            groFile.write(line)
        
        box = max([max(self.coordsX),max(self.coordsY),max(self.coordsZ)])
        boxX=boxY=boxZ=box+(box*0.2) 
        text = "%11.5f %11.5f %11.5f\n" % (boxX, boxY, boxZ )
        

        groFile.write(text)
    
        groFile.close()
    


    def writeTop(self):

        topText = []
        molName = self.name

        headTop = """
# include "martini_v3.0.0.itp"
# include "martini_v3.0.0_solvents_v1.itp"
"""

        headMoleculetype = """
[ moleculetype ]
;name            nrexcl
%-16s 1
"""

        headAtomTypes = """
[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
"""

        headBonds = """
[ bonds ]
;   ai     aj funct   r             k
"""

        headAngles = """
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""
        headProDih = """
[ dihedrals ] 
;    i      j      k      l   func    C0         C1         C2         C3         C4         C5
"""
        headImp = """
; impropers  
"""

#        topText.append(headTop)
        topText.append(headMoleculetype % molName)
        topText.append(headAtomTypes)
        for row in self.atoms:
            line = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot %1.3f\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])
            topText.append(line)

        topText.append(headBonds)
        for row in self.bonds:
            line = "%6i %6i %3i %13.4e %13.4e\n" % (row[0], row[1], 1, row[2], row[3])
            topText.append(line)

        topText.append(headAngles)
        for row in self.angles:
            line = "%6i %6i %6i %6i %13.4e %13.4e\n" % (row[0], row[1], row[2], row[3], row[4], row[5])
            topText.append(line)

        topText.append(headProDih)
        for row in self.dihedrals:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
            topText.append(line)

        topText.append(headImp)
        for row in self.impropers:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
            topText.append(line)

        footerTop = """
[ system ]
; Name
UNNAMED

[ molecules ]
; Compound        #mols
%-16s             1
"""

#        topText.append(footerTop % molName)

        fileName = molName + ".itp"
        topFile = open(fileName, "w")
        topFile.writelines(topText)
        topFile.close()

    
def GLCA14(units_num):

    name = "GLCA14"
    at = 4
    

    atom = []
    atom.append([1, "TP3", "AMY", "T1", 0, 36])
    atom.append([2, "TN4", "AMY", "T2", 0, 36])
    atom.append([3, "P3", "AMY", "R3", 0, 72])
    atom.append([4, "SN4", "AMY", "S4", 0, 54])
    
    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])

                
    bondIn = []
    bondIn.append([1, 4, 0.251, 18000])
    bondIn.append([2, 3, 0.280, 32000])
    bondIn.append([2, 4, 0.292, 60000])
    bondIn.append([3, 4, 0.279, 34000])

    bondConnect = []
    bondConnect.append([2, 8, 0.280, 32000])

    angleIn = []
    angleIn.append([1, 4, 2, 2,  86, 340])
    angleIn.append([1, 4, 3, 2, 142, 580])
    
    angleConnect = []
    angleConnect.append([3, 2, 8,  2, 103, 310]) 
    angleConnect.append([4, 2, 8, 10, 103, 180]) 
    angleConnect.append([1, 2, 8,  2, 105, 210]) 
    angleConnect.append([2, 8, 5,  2, 107,  30]) 
    angleConnect.append([2, 8, 6, 10, 148, 240]) 
    angleConnect.append([2, 8, 7,  2,  99, 120]) 
    
    dihIn = []
    
    
    dihConnect = []
    dihConnect.append([4, 2, 8, 6, 1, 161, -12.8, 1])
    dihConnect.append([3, 2, 8, 7, 1, 165, -2.2, 1])
    dihConnect.append([1, 3, 7, 5, 1, -72, -3.6, 5])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 9, 200])
    
    improperConnect = []
    improperConnect.append([8, 2, 7, 6, 2, 22.3, 170])
    improperConnect.append([2, 8, 4, 3, 2, -67, 300])



    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    axis = [1,0,0]
    for i in range(units_num):
        if (i+1) % 2 == 1:
            theta = 0
        elif (i+1) % 2 == 0:
            theta = 3.1415926536 # 180
        for row in coord:
            #row = rotation_matrix(axis, theta) @ row
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
                mol.gets_impropers(0, row)

    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()
        

def GLCB12(units_num):

    name = "GLCB12"
    at = 4
    dist = 0.55

    atom = []
    atom.append([1, "TP3", "GLC2", "T1", 0, 36])
    atom.append([2, "TN4", "GLC2", "T2", 0, 36])
    atom.append([3, "P3",  "GLC2", "R3", 0, 72])
    atom.append([4, "SN4", "GLC2", "S4", 0, 54])

    coord = []

    coord.append([0.454, 0.200, 0.000])
    coord.append([0.275, 0.401, 0.000])
    coord.append([0.000, 0.400, 0.000])
    coord.append([0.202, 0.215, 0.000])

    bondIn = []
    bondIn.append([1, 4, 0.265, 10000])
    bondIn.append([2, 3, 0.263, 50000])
    bondIn.append([2, 4, 0.290, 50000])
    bondIn.append([3, 4, 0.272, 24000])

    bondConnect = []
    bondConnect.append([2, 7, 0.274, 6000])
    bondConnect.append([4, 7, 0.530, 600])


    angleIn = []
    angleIn.append([1, 4, 2, 2, 80, 200])
    angleIn.append([1, 4, 3, 2, 132, 400])
    

    angleConnect = []
    angleConnect.append([3, 2, 7, 2, 88, 30])
    angleConnect.append([1, 2, 7, 10, 165, 60])
    angleConnect.append([2, 7, 6, 2, 97, 50])
    angleConnect.append([2, 7, 8, 10, 148, 40]) 


    dihIn = []

    dihConnect = []
    dihConnect.append([3, 2, 7, 6, 1, -141, -22, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 11, 150])

    improperConnect = []
    improperConnect.append([7, 6, 8, 2, 2, 15, 300])
    improperConnect.append([2, 4, 3, 7, 2, 1, 200])

    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54
    
    for i in range(units_num):
        for row in coord:
            mol.coordsX.append(row[0]+i*dist)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()
        


def GLCA16(units_num):

    name = "GLCA16"
    at = 4 
    

    atom = []
    atom.append([1, "TP3", "GLC6", "T1", 0, 36])
    atom.append([2, "TN4", "GLC6", "T2", 0, 36])
    atom.append([3, "P3", "GLC6", "R3", 0, 72])
    atom.append([4, "SN4", "GLC6", "S4", 0, 54])

    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])


    bondIn = []
    bondIn.append([1, 4, 0.251, 21000])
    bondIn.append([2, 3, 0.280, 35000])
    bondIn.append([2, 4, 0.322, 60000])
    bondIn.append([3, 4, 0.289, 24000])

    bondConnect = []
    bondConnect.append([2, 5, 0.225, 13000])
    bondConnect.append([3, 7, 0.82, 2400])

    bondConnect3 = []
    bondConnect3.append([1, 9, 0.79, 180])
    
    angleIn = []
    angleIn.append([1, 4, 2, 2, 80, 600])
    angleIn.append([1, 4, 3, 2, 135, 560])

    angleConnect = []
    angleConnect.append([3, 2, 5, 2, 109, 180])
    angleConnect.append([4, 2, 5, 2, 105, 100])
    angleConnect.append([2, 5, 6, 10, 122, 28])

       

    dihIn = []
    
    dihConnect = []
    dihConnect.append([3, 2, 5, 6, 1, 155, -3, 2])
    dihConnect.append([2, 5, 6, 7, 1, 170, -3, 5])
    dihConnect.append([2, 5, 6, 7, 1, 60, 3, 1])
    dihConnect.append([2, 4, 8, 6, 1, -25, 5, 1])
    dihConnect.append([2, 4, 8, 6, 1, 165, -4, 3])

    improperIn = []
    improperIn.append([4, 2, 3, 1, 2, -7, 250])

    improperConnect = []
    improperConnect.append([2, 3, 5, 4, 2, -74, 110])


    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    axis = [1,0,0]
    for i in range(units_num):
        for row in coord:
            theta = pickAngle(i)
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    if bondConnect3:
        for row in bondConnect3:
            for i in range(units_num - 2):
                mol.gets_bonds(i, row)


    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()



def GLCB13(units_num):

    name = "GLCB13"
    at = 4
    dist = 0.55

    atom = []
    atom.append([1, "TP3", "CUR", "T1", 0, 36])
    atom.append([2, "TN4", "CUR", "T2", 0, 36])
    atom.append([3, "P3", "CUR", "R3", 0, 72])
    atom.append([4, "SN4", "CUR", "S4", 0, 54])

    coord = []

    coord.append([0.454, 0.200, 0.000])
    coord.append([0.275, 0.401, 0.000])
    coord.append([0.000, 0.400, 0.000])
    coord.append([0.202, 0.215, 0.000])

    bondIn = []
    bondIn.append([1, 4, 0.268, 10000])
    bondIn.append([2, 3, 0.240, 38000])
    bondIn.append([2, 4, 0.287, 50000])
    bondIn.append([3, 4, 0.287, 24000])

    bondConnect = []
    bondConnect.append([2, 7, 0.250, 12000])


    angleIn = []
    angleIn.append([1, 4, 2, 2, 78, 220])
    angleIn.append([1, 4, 3, 2, 123, 400])

    angleConnect = []
    angleConnect.append([3, 2, 7, 2, 75, 60]) 
    angleConnect.append([1, 2, 7, 10, 148, 80]) 
    angleConnect.append([2, 7, 6, 2, 128, 110]) 
    angleConnect.append([2, 7, 8, 2, 65, 40]) 

    dihIn = []

    dihConnect = []
    dihConnect.append([3, 2, 7, 6, 1, -174, -20, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 15, 120])

    improperConnect = []
    improperConnect.append([7, 2, 8, 6, 2, 22, 200])
    improperConnect.append([2, 7, 3, 4, 2, 3, 200])

    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    
    for i in range(units_num):
        for row in coord:
            mol.coordsX.append(row[0]+i*dist)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 
    

    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()


def GLCB14(units_num):

    name = "GLCB14"
    at = 4
    

    atom = []
    atom.append([1, "TP3", "CELL", "T1", 0, 36])
    atom.append([2, "TN4", "CELL", "T2", 0, 36])
    atom.append([3, "P3", "CELL", "R3", 0, 72])
    atom.append([4, "SN4", "CELL", "S4", 0, 54])
    
    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])

    

    bondIn = []
    bondIn.append([1, 4, 0.250, 14100])
    bondIn.append([2, 3, 0.2680, 37500])
    bondIn.append([2, 4, 0.2570, 53200])
    bondIn.append([3, 4, 0.2730, 27000])

    bondConnect = []
    bondConnect.append([2, 8, 0.267, 7500])
    bondConnect.append([2, 6, 0.520, 16300])
    bondConnect.append([4, 8, 0.542, 3770])

    angleIn = []
    angleIn.append([1, 4, 2, 2, 91, 220])
    angleIn.append([1, 4, 3, 10, 143, 159])
    
    angleConnect = []
    angleConnect.append([3, 2, 8, 2, 115, 245]) 
    angleConnect.append([1, 2, 8, 2, 127, 350]) 
    angleConnect.append([2, 8, 5, 2, 123, 16]) 
    angleConnect.append([2, 8, 7, 2, 93, 52]) 
    
    dihIn = []
    

    dihConnect = []
    dihConnect.append([3, 2, 8, 7, 1, -135 ,-35, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 9, 200])

    improperConnect = []
    improperConnect.append([2, 3, 8, 4, 2, 9, 229])
    improperConnect.append([8, 2, 7, 6, 2, 11, 212])   


    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54
        

    axis = [1,0,0]
    for i in range(units_num):
        if (i+1) % 2 == 1:
            theta = 0
        elif (i+1) % 2 == 0:
            theta = 3.1415926536 # 180
        for row in coord:
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
                mol.gets_impropers(0, row)

    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()



def HA(units_num):
    name = "HA"
    at = 9
    dist = 1.14

    atom = []
    atom.append([1, "TP3r ", "HA", "T1", 0, 31])
    atom.append([2, "TN4r ", "HA", "T2", 0, 29])
    atom.append([3, " P2r ", "HA", "R3", 0, 57])
    atom.append([4, "SN4r ", "HA", "S4", 0, 43])
    atom.append([5, "SN5ar", "HA", "S5", 0, 43])
    atom.append([6, "SQ5n ", "HA", "S6", -1, 44])
    atom.append([7, "TN4r ", "HA", "T7", 0, 29])
    atom.append([8, " P3r ", "HA", "R8", 0, 60])
    atom.append([9, "SN4r ", "HA", "S9", 0, 42])

    coord = []

    coord.append([1.609, 3.524, 3.389])
    coord.append([1.789, 3.762, 3.565])
    coord.append([1.655, 3.972, 3.616])
    coord.append([1.503, 3.743, 3.528])
    coord.append([1.789, 4.238, 3.723])
    coord.append([2.088, 4.070, 3.636])
    coord.append([2.311, 3.789, 3.533])
    coord.append([2.116, 3.597, 3.514])
    coord.append([2.062, 3.859, 3.550])

    bondIn = []
    bondIn.append([1, 4, 0.274, 11700])
    bondIn.append([6, 9, 0.233, 65000])
    bondIn.append([2, 3, 0.232, 73000])
    bondIn.append([7, 8, 0.261, 26000])
    bondIn.append([2, 4, 0.289, 36000])
    bondIn.append([7, 9, 0.257, 50000])
    bondIn.append([3, 4, 0.284, 42000])
    bondIn.append([8, 9, 0.278, 28000])
    bondIn.append([3, 5, 0.330, 38000])
    bondIn.append([2, 9, 0.255, 20000])
    bondIn.append([4, 5, 0.615, 3000])
    bondIn.append([4, 7, 0.799, 10000])

    bondConnect = []
    bondConnect.append([7, 12, 0.292, 400])
    bondConnect.append([7, 11, 0.473, 5000])
    bondConnect.append([9, 13, 0.573, 4000])
    bondConnect.append([2, 12, 0.802, 5000])

    angleIn = []
    angleIn.append([1, 4, 2, 2, 78, 345])
    angleIn.append([6, 9, 7, 2, 93, 630])
    angleIn.append([1, 4, 3, 2, 126, 270])
    angleIn.append([6, 9, 8, 10, 146, 345])
    angleIn.append([5, 3, 2, 2, 126, 500])
    angleIn.append([3, 2, 9, 2, 93, 150])
    angleIn.append([2, 9, 8, 2, 116, 150])
    angleIn.append([2, 9, 6, 2, 80, 150])

    angleConnect = []
    angleConnect.append([8, 7, 12, 2, 110, 250])
    angleConnect.append([7, 12, 11, 2, 120, 200])
    angleConnect.append([7, 12, 13, 2, 66, 150])

    dihIn = []
    dihIn.append([1, 3, 8, 6, 1, -138, -32, 1])

    dihConnect = []
    dihConnect.append([8, 7, 12, 11, 1, -165, -60, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 16, 150])
    improperIn.append([5, 3, 2, 4, 2, -173, 170])
    improperIn.append([9, 8, 7, 6, 2, 13, 350])
    improperIn.append([2, 3, 9, 4, 2, 5, 170])
    improperIn.append([9, 2, 8, 7, 2, 13, 300])

    improperConnect = []
    improperConnect.append([7, 8, 9, 12, 2, 1, 500])
    improperConnect.append([12, 11, 13, 7, 2, -18, 400])

    mol = Polysaccharide(name, at)

    for i in range(units_num):
        for row in atom:
            mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4r "
    mol.atoms[-3][4] = "S7"
    mol.atoms[-3][7] = 46
    mol.atoms[2][7] = 58

    for i in range(units_num):
        for row in coord:
            mol.coordsX.append(row[0] + i * dist + 5) #added + 5
            mol.coordsY.append(row[1] + 15) #added + 15
            mol.coordsZ.append(row[2] + 15) #added + 15

    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num - 1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num - 1):
                mol.gets_dihedrals(i, row)

    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num - 1):
                mol.gets_impropers(i, row)

    mol.writeGro()
    mol.writeTop()



def _vector_angle_deg(a, b, c):
    ba = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    bc = np.asarray(c, dtype=float) - np.asarray(b, dtype=float)
    nba = np.linalg.norm(ba)
    nbc = np.linalg.norm(bc)
    if nba == 0 or nbc == 0:
        return 180.0
    cosang = np.dot(ba, bc) / (nba * nbc)
    cosang = np.clip(cosang, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosang)))


def _dihedral_deg(p1, p2, p3, p4):
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    p3 = np.asarray(p3, dtype=float)
    p4 = np.asarray(p4, dtype=float)

    b0 = p2 - p1
    b1 = p3 - p2
    b2 = p4 - p3

    nb1 = np.linalg.norm(b1)
    if nb1 == 0:
        return 180.0
    b1 = b1 / nb1

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    nv = np.linalg.norm(v)
    nw = np.linalg.norm(w)
    if nv == 0 or nw == 0:
        return 180.0

    v /= nv
    w /= nw

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def _rotation_from_vectors(vec1, vec2):
    v1 = np.asarray(vec1, dtype=float)
    v2 = np.asarray(vec2, dtype=float)

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return np.eye(3)

    v1 /= n1
    v2 /= n2
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)

    if np.linalg.norm(cross) < 1e-12:
        if dot > 0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(v1[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0])
        axis = np.cross(v1, axis)
        axis /= np.linalg.norm(axis)
        return rotation_matrix(axis, np.pi)

    axis = cross / np.linalg.norm(cross)
    angle = float(np.arccos(np.clip(dot, -1.0, 1.0)))
    return rotation_matrix(axis, angle)


def _ha_osa_templates():
    ha_atoms = [
        [1, "TP3r ", "HA", "T1", 0, 31],
        [2, "TN4r ", "HA", "T2", 0, 29],
        [3, " P2r ", "HA", "R3", 0, 57],
        [4, "SN4r ", "HA", "S4", 0, 43],
        [5, "SN5ar", "HA", "S5", 0, 43],
        [6, "SQ5n ", "HA", "S6", -1, 44],
        [7, "TN4r ", "HA", "T7", 0, 29],
        [8, " P3r ", "HA", "R8", 0, 60],
        [9, "SN4r ", "HA", "S9", 0, 42],
    ]

    ha_coords = {
        1: np.array([1.609, 3.524, 3.389], dtype=float),
        2: np.array([1.789, 3.762, 3.565], dtype=float),
        3: np.array([1.655, 3.972, 3.616], dtype=float),
        4: np.array([1.503, 3.743, 3.528], dtype=float),
        5: np.array([1.789, 4.238, 3.723], dtype=float),
        6: np.array([2.088, 4.070, 3.636], dtype=float),
        7: np.array([2.311, 3.789, 3.533], dtype=float),
        8: np.array([2.116, 3.597, 3.514], dtype=float),
        9: np.array([2.062, 3.859, 3.550], dtype=float),
    }

    ha_bond_in = [
        [1, 4, 0.274, 11700],
        [6, 9, 0.233, 65000],
        [2, 3, 0.232, 73000],
        [7, 8, 0.261, 26000],
        [2, 4, 0.289, 36000],
        [7, 9, 0.257, 50000],
        [3, 4, 0.284, 42000],
        [8, 9, 0.278, 28000],
        [3, 5, 0.330, 38000],
        [2, 9, 0.255, 20000],
        [4, 5, 0.615, 3000],
        [4, 7, 0.799, 10000],
    ]

    ha_bond_connect = [
        [7, 12, 0.292, 400],
        [7, 11, 0.473, 5000],
        [9, 13, 0.573, 4000],
        [2, 12, 0.802, 5000],
    ]

    ha_angle_in = [
        [1, 4, 2, 2, 78, 345],
        [6, 9, 7, 2, 93, 630],
        [1, 4, 3, 2, 126, 270],
        [6, 9, 8, 10, 146, 345],
        [5, 3, 2, 2, 126, 500],
        [3, 2, 9, 2, 93, 150],
        [2, 9, 8, 2, 116, 150],
        [2, 9, 6, 2, 80, 150],
    ]

    ha_angle_connect = [
        [8, 7, 12, 2, 110, 250],
        [7, 12, 11, 2, 120, 200],
        [7, 12, 13, 2, 66, 150],
    ]

    ha_dih_in = [
        [1, 3, 8, 6, 1, -138, -32, 1],
    ]

    ha_dih_connect = [
        [8, 7, 12, 11, 1, -165, -60, 1],
    ]

    ha_improper_in = [
        [4, 3, 2, 1, 2, 16, 150],
        [5, 3, 2, 4, 2, -173, 170],
        [9, 8, 7, 6, 2, 13, 350],
        [2, 3, 9, 4, 2, 5, 170],
        [9, 2, 8, 7, 2, 13, 300],
    ]

    ha_improper_connect = [
        [7, 8, 9, 12, 2, 1, 500],
        [12, 11, 13, 7, 2, -18, 400],
    ]

    osa_atoms = [
        [1, "SN4a ", "OSA", "E1", 0.0, 58.0],
        [2, "TC3  ", "OSA", "E2", 0.0, 14.0],
        [3, "TP2a ", "OSA", "E3", -1.0, 45.0],
        [4, "SC4  ", "OSA", "E4", 0.0, 56.0],
        [5, "C1   ", "OSA", "E5", 0.0, 58.0],
    ]

    osa_coords = {
        1: np.array([0.04500000, -0.00200000, 0.01800000], dtype=float),
        2: np.array([0.31500000, -0.00200000, 0.01800000], dtype=float),
        3: np.array([0.47500000,  0.27512813, 0.01800000], dtype=float),
        4: np.array([0.45500000, -0.13897230, 0.21809645], dtype=float),
        5: np.array([0.63915265, -0.24626199, 0.42922845], dtype=float),
    }

    osa_bonds = [
        [1, 2, 0.2700, 25000.0],
        [2, 3, 0.3200, 22000.0],
        [2, 4, 0.2800, 20000.0],
        [4, 5, 0.3000, 12000.0],
    ]

    osa_angles = [
        [1, 2, 3, 2, 120.0, 80.0],
        [1, 2, 4, 2, 120.0, 80.0],
        [3, 2, 4, 2, 100.0, 70.0],
        [2, 4, 5, 2, 170.0, 50.0],
    ]

    return {
        "ha_atoms": ha_atoms,
        "ha_coords": ha_coords,
        "ha_bond_in": ha_bond_in,
        "ha_bond_connect": ha_bond_connect,
        "ha_angle_in": ha_angle_in,
        "ha_angle_connect": ha_angle_connect,
        "ha_dih_in": ha_dih_in,
        "ha_dih_connect": ha_dih_connect,
        "ha_improper_in": ha_improper_in,
        "ha_improper_connect": ha_improper_connect,
        "osa_atoms": osa_atoms,
        "osa_coords": osa_coords,
        "osa_bonds": osa_bonds,
        "osa_angles": osa_angles,
        "ha_dist": 1.14,
        "ha_shift": np.array([5.0, 15.0, 15.0], dtype=float),
    }


def _ha_absolute_coord(repeat_index, local_id, templates):
    base = templates["ha_coords"][local_id]
    shift = templates["ha_shift"] + np.array([(repeat_index - 1) * templates["ha_dist"], 0.0, 0.0], dtype=float)
    return base + shift


def _resolve_repeat_local(base_repeat, raw_index, unit_at):
    repeat_index = base_repeat + (raw_index - 1) // unit_at
    local_id = ((raw_index - 1) % unit_at) + 1
    return repeat_index, local_id


def _place_osa_fragment(repeat_index, templates):
    osa_local = templates["osa_coords"]
    e1_ref = osa_local[1]
    centered = {idx: coord - e1_ref for idx, coord in osa_local.items()}

    t1_pos = _ha_absolute_coord(repeat_index, 1, templates)
    s4_pos = _ha_absolute_coord(repeat_index, 4, templates)

    outward = t1_pos - s4_pos
    source = centered[2] - centered[1]

    rot = _rotation_from_vectors(source, outward)
    aligned = {idx: np.dot(rot, coord) for idx, coord in centered.items()}

    local_neighborhood = []
    for rep in range(max(1, repeat_index - 1), repeat_index + 2):
        for local_id in range(1, 10):
            if rep == repeat_index and local_id == 1:
                continue
            local_neighborhood.append(_ha_absolute_coord(rep, local_id, templates))

    axis = outward / np.linalg.norm(outward)
    best_score = -1.0e30
    best_coords = None

    for phi in np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False):
        spin = rotation_matrix(axis, phi)
        trial = {idx: np.dot(spin, coord) + t1_pos for idx, coord in aligned.items()}

        distances = [
            np.linalg.norm(trial[idx] - neigh)
            for idx in [2, 3, 4, 5]
            for neigh in local_neighborhood
        ]
        min_dist = min(distances)

        penalty = 0.0
        for d in distances:
            if d < 0.60:
                penalty += (0.60 - d) ** 2 * 200.0

        score = min_dist + 0.20 * float(np.mean(distances)) - penalty

        if score > best_score:
            best_score = score
            best_coords = trial

    return best_coords


def OSA_HA(units_num, ds_percent, seed):
    templates = _ha_osa_templates()
    substitution_count = math.floor((ds_percent / 100.0) * units_num)

    if substitution_count <= 0:
        print("Zero substitutions selected from the requested DS. Generating plain HA.")
        HA(units_num)
        return

    rng = random.Random(seed)
    substituted_repeats = sorted(rng.sample(range(1, units_num + 1), substitution_count))
    substituted_set = set(substituted_repeats)

    mol = Polysaccharide("OSA_HA", 0)
    qtot = 0.0
    next_atom_index = 1

    ha_map = {}
    atom_coords = {}

    for repeat_index in range(1, units_num + 1):
        for row in templates["ha_atoms"]:
            local_id = row[0]
            if repeat_index in substituted_set and local_id == 1:
                continue

            bead_type = row[1]
            bead_name = row[3]
            mass = row[5]
            if repeat_index == 1 and local_id == 3:
                mass = 58
            if repeat_index == units_num and local_id == 7:
                bead_type = "SN4r "
                bead_name = "S7"
                mass = 46

            qtot += row[4]
            mol.atoms.append([next_atom_index, bead_type, repeat_index, row[2], bead_name, next_atom_index, row[4], mass, qtot])
            coord = _ha_absolute_coord(repeat_index, local_id, templates)
            mol.coordsX.append(coord[0])
            mol.coordsY.append(coord[1])
            mol.coordsZ.append(coord[2])

            ha_map[(repeat_index, local_id)] = next_atom_index
            atom_coords[next_atom_index] = coord
            next_atom_index += 1

    osa_map = {}
    osa_residue_index = units_num + 1

    for repeat_index in substituted_repeats:
        placed = _place_osa_fragment(repeat_index, templates)
        osa_map[repeat_index] = {}
        for row in templates["osa_atoms"]:
            local_id = row[0]
            bead_type = row[1]
            bead_name = row[3]
            charge = row[4]
            mass = row[5]

            qtot += charge
            mol.atoms.append([next_atom_index, bead_type, osa_residue_index, row[2], bead_name, next_atom_index, charge, mass, qtot])
            coord = placed[local_id]
            mol.coordsX.append(coord[0])
            mol.coordsY.append(coord[1])
            mol.coordsZ.append(coord[2])

            osa_map[repeat_index][local_id] = next_atom_index
            atom_coords[next_atom_index] = coord
            next_atom_index += 1

        osa_residue_index += 1

    def map_atom(repeat_index, local_id):
        if repeat_index in substituted_set and local_id == 1:
            return osa_map[repeat_index][1]
        return ha_map.get((repeat_index, local_id))

    def add_bond(ai, aj, r, k):
        mol.bonds.append([ai, aj, r, k])

    def add_angle(ai, aj, ak, funct, theta, cth):
        mol.angles.append([ai, aj, ak, funct, theta, cth])

    def add_dihedral(ai, aj, ak, al, funct, c0, c1, c2):
        mol.dihedrals.append([ai, aj, ak, al, funct, c0, c1, c2])

    def add_improper(ai, aj, ak, al, funct, angle, force):
        mol.impropers.append([ai, aj, ak, al, funct, angle, force])

    for row in templates["ha_bond_in"]:
        for repeat_index in range(1, units_num + 1):
            ai = map_atom(repeat_index, row[0])
            aj = map_atom(repeat_index, row[1])
            if ai is not None and aj is not None:
                add_bond(ai, aj, row[2], row[3])

    for row in templates["ha_bond_connect"]:
        for repeat_index in range(1, units_num):
            rep_i, local_i = _resolve_repeat_local(repeat_index, row[0], 9)
            rep_j, local_j = _resolve_repeat_local(repeat_index, row[1], 9)
            ai = map_atom(rep_i, local_i)
            aj = map_atom(rep_j, local_j)
            if ai is not None and aj is not None:
                add_bond(ai, aj, row[2], row[3])

    for row in templates["ha_angle_in"]:
        for repeat_index in range(1, units_num + 1):
            ai = map_atom(repeat_index, row[0])
            aj = map_atom(repeat_index, row[1])
            ak = map_atom(repeat_index, row[2])
            if ai is not None and aj is not None and ak is not None:
                add_angle(ai, aj, ak, row[3], row[4], row[5])

    for row in templates["ha_angle_connect"]:
        for repeat_index in range(1, units_num):
            rep_i, local_i = _resolve_repeat_local(repeat_index, row[0], 9)
            rep_j, local_j = _resolve_repeat_local(repeat_index, row[1], 9)
            rep_k, local_k = _resolve_repeat_local(repeat_index, row[2], 9)
            ai = map_atom(rep_i, local_i)
            aj = map_atom(rep_j, local_j)
            ak = map_atom(rep_k, local_k)
            if ai is not None and aj is not None and ak is not None:
                add_angle(ai, aj, ak, row[3], row[4], row[5])

    for row in templates["ha_dih_in"]:
        for repeat_index in range(1, units_num + 1):
            ai = map_atom(repeat_index, row[0])
            aj = map_atom(repeat_index, row[1])
            ak = map_atom(repeat_index, row[2])
            al = map_atom(repeat_index, row[3])
            if ai is not None and aj is not None and ak is not None and al is not None:
                add_dihedral(ai, aj, ak, al, row[4], row[5], row[6], row[7])

    for row in templates["ha_dih_connect"]:
        for repeat_index in range(1, units_num):
            rep_i, local_i = _resolve_repeat_local(repeat_index, row[0], 9)
            rep_j, local_j = _resolve_repeat_local(repeat_index, row[1], 9)
            rep_k, local_k = _resolve_repeat_local(repeat_index, row[2], 9)
            rep_l, local_l = _resolve_repeat_local(repeat_index, row[3], 9)
            ai = map_atom(rep_i, local_i)
            aj = map_atom(rep_j, local_j)
            ak = map_atom(rep_k, local_k)
            al = map_atom(rep_l, local_l)
            if ai is not None and aj is not None and ak is not None and al is not None:
                add_dihedral(ai, aj, ak, al, row[4], row[5], row[6], row[7])

    for row in templates["ha_improper_in"]:
        for repeat_index in range(1, units_num + 1):
            ai = map_atom(repeat_index, row[0])
            aj = map_atom(repeat_index, row[1])
            ak = map_atom(repeat_index, row[2])
            al = map_atom(repeat_index, row[3])
            if ai is not None and aj is not None and ak is not None and al is not None:
                add_improper(ai, aj, ak, al, row[4], row[5], row[6])

    for row in templates["ha_improper_connect"]:
        for repeat_index in range(1, units_num):
            rep_i, local_i = _resolve_repeat_local(repeat_index, row[0], 9)
            rep_j, local_j = _resolve_repeat_local(repeat_index, row[1], 9)
            rep_k, local_k = _resolve_repeat_local(repeat_index, row[2], 9)
            rep_l, local_l = _resolve_repeat_local(repeat_index, row[3], 9)
            ai = map_atom(rep_i, local_i)
            aj = map_atom(rep_j, local_j)
            ak = map_atom(rep_k, local_k)
            al = map_atom(rep_l, local_l)
            if ai is not None and aj is not None and ak is not None and al is not None:
                add_improper(ai, aj, ak, al, row[4], row[5], row[6])

    for repeat_index in substituted_repeats:
        e = osa_map[repeat_index]
        s4 = ha_map[(repeat_index, 4)]

        for row in templates["osa_bonds"]:
            add_bond(e[row[0]], e[row[1]], row[2], row[3])

        for row in templates["osa_angles"]:
            add_angle(e[row[0]], e[row[1]], e[row[2]], row[3], row[4], row[5])

        # No extra graft angle here
        # No graft cross-dihedral here

    mol.writeGro()
    mol.writeTop()

    print("Generated OSA_HA with the following summary:")
    print("Total HA disaccharide repeats: %d" % units_num)
    print("Requested DS (percent): %s" % ds_percent)
    print("Applied substitutions: %d" % substitution_count)
    print("Substituted repeat indices (1-based): %s" % substituted_repeats)
    print("Total atoms: %d" % len(mol.atoms))
    print("Total bonds: %d" % len(mol.bonds))
    print("Total angles: %d" % len(mol.angles))
    print("Total dihedrals: %d" % len(mol.dihedrals))
    print("Total impropers: %d" % len(mol.impropers))


def main():

    units_name = ["alpha(1-4)", "beta(1-2)", "alpha(1-6)", "beta(1-4)", "beta(1-3)", "hyaluronic_acid(HA),monomer=dimer[GlcNAc–beta(1-4)–beta-GlcA]", "OSA_hyaluronic_acid(OSA_HA)"]
    print("Select linkage type:")
    for i, name in enumerate(units_name, 1):
        print("%d %s" % (i, name))
    sugar = int(input())
    print("You chose %s" % (units_name[sugar - 1]))
    while True:
        print("Insert number of residue in a homopolysaccharide:")
        units_nr = int(input())
        if(units_nr==1):
            print("The number of residues must be larger than 1")
            continue
        else:
            break
    if(sugar == 1):
        GLCA14(units_nr)
    if(sugar == 2):
        GLCB12(units_nr)
    if(sugar == 3):
        GLCA16(units_nr)
    if(sugar == 4):
        GLCB14(units_nr)
    if(sugar == 5):
        GLCB13(units_nr)
    if (sugar == 6):
        HA(units_nr)
    if (sugar == 7):
        print("Insert DS percentage for OSA substitution per HA disaccharide:")
        ds_percent = float(input())
        print("Insert random seed:")
        seed = int(input())
        OSA_HA(units_nr, ds_percent, seed)

if __name__ == "__main__":
    main()
