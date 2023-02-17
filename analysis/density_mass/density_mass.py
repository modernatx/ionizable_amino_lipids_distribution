import numpy as np
import math  as m
import MDAnalysis as mda
import fnmatch,os
from MDAnalysis import transformations

def trajectory_multiple(dcdfile, start, end):
    list_traj = []
    for i in range(start, end):
        list_traj.append(dcdfile + "." + str(i) + ".dcd")
    return list_traj

LIST1 = ["10/1/equ", "20/1/equ_1", "30/1/equ", "40/1/equ_1", "50/1/equ_1"]
LIST2 = ["10/1/"   , "20/1/"     , "30/1/"   , "40/1/"     , "50/1/"]
STEP  = 1

for p in range(0, 5): 
    
    psffile   = "/home/sdehghanighahnaviyeh/SM86_systems/" + str(LIST2[p]) + "/output/ionized.psf"
    dcdfile   = "/home/sdehghanighahnaviyeh/SM86_systems/" + str(LIST1[p]) + "/reduced.dcd"

    u = mda.Universe(psffile, dcdfile)
    N_frames = len(u.trajectory)
    X_length = u.dimensions[0]
    Y_length = u.dimensions[1]
    START = int(N_frames * 0.666666666666)

    print ("trajectory is read :-)")
    print ("Number of frames: " + str(N_frames))
    print ("The dimension of the box in x-y plane: x = " + str(X_length) + " and y = " + str(Y_length))
    
    OUTPUT1 = open("density-" + str(p) + ".txt" , "w")
    
    z1 = np.array([])
    z2 = np.array([])
    z3 = np.array([])
    z4 = np.array([])
    z5 = np.array([])

    for ts in u.trajectory[START::STEP]:
        if ts.frame % 1 == 0:
            print (str(ts.frame) + " out of " + str(N_frames))
        phosphates = u.select_atoms("name P31")
        COM = phosphates.center_of_geometry()
        sel1 = u.select_atoms("name N and resname RP SP")
        R1 = sel1.positions
        for i in range(0, len(R1)):
            z1    = (np.append(z1, R1[i][2] - COM))
        sel2 = u.select_atoms("name N and resname RN SN")
        R2 = sel2.positions
        for i in range(0, len(R2)):
            z2    = (np.append(z2, R2[i][2] - COM))
        sel3 = u.select_atoms("name P31")
        R3 = sel3.positions
        for i in range(0, len(R3)):
            z3    = (np.append(z3, R3[i][2] - COM))
        sel4 = u.select_atoms("name N")
        R4 = sel4.positions
        for i in range(0, len(R4)):
            z4    = (np.append(z4, R4[i][2] - COM))

    hist1, bin_edges = np.histogram(z1, bins=100, range=(-70, 70))
    hist2, bin_edges = np.histogram(z2, bins=100, range=(-70, 70))
    hist3, bin_edges = np.histogram(z3, bins=100, range=(-70, 70))
    hist4, bin_edges = np.histogram(z4, bins=100, range=(-70, 70))

    for i in range(0, len(bin_edges) - 1):
        OUTPUT1.write('{}     '.format(bin_edges[i]))
        OUTPUT1.write('{}     '.format(hist1[i]))
        OUTPUT1.write('{}     '.format(hist2[i]))
        OUTPUT1.write('{}     '.format(hist3[i]))
        OUTPUT1.write('{}   \n'.format(hist4[i]))
