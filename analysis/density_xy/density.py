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
    
    OUTPUT1 = open("density-u" + str(p) + ".txt" , "w")
    OUTPUT2 = open("density-l" + str(p) + ".txt" , "w")

    for ts in u.trajectory[START::STEP]:
        if ts.frame % 1 == 0:
            print (str(ts.frame) + " out of " + str(N_frames))

        sel1 = u.select_atoms("name N and resname RP SP and prop z > 0")
        R1 = sel1.positions
        
        for i in range(0, len(R1)):
            OUTPUT1.write('{}     '.format(R1[i][0]))
            OUTPUT1.write('{}     '.format(R1[i][1]))
            OUTPUT1.write('{}   \n'.format(R1[i][2]))
            
        sel2 = u.select_atoms("name N and resname RP SP and prop z < 0")
        R2 = sel2.positions
        
        for i in range(0, len(R2)):
            OUTPUT2.write('{}     '.format(R2[i][0]))
            OUTPUT2.write('{}     '.format(R2[i][1]))
            OUTPUT2.write('{}   \n'.format(R2[i][2]))
