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
#     num_files = len(fnmatch.filter(os.listdir("/home/sdehghanighahnaviyeh/SM86_systems/" + str(LIST1[p])), 'equ.*.dcd'))
#     list_traj = trajectory_multiple(dcdfile, 0, num_files)
#     print (list_traj)
    u = mda.Universe(psffile, dcdfile)
    N_frames = len(u.trajectory)
    X_length = u.dimensions[0]
    Y_length = u.dimensions[1]
    print ("trajectory is read :-)")
    print ("Number of frames: " + str(N_frames))
    print ("The dimension of the box in x-y plane: x = " + str(X_length) + " and y = " + str(Y_length))
#     ag = u.atoms
#     transform = mda.transformations.unwrap(ag)
#     u.trajectory.add_transformations(transform)

    OUTPUT1 = open("coor-n-" + str(p) + ".txt" , "w")
    OUTPUT2 = open("coor-p-" + str(p) + ".txt" , "w")
    START = int(N_frames * 0.666666666666)

    for ts in u.trajectory[START::STEP]:
        if ts.frame % 10 == 0:
            print (str(ts.frame) + " out of " + str(N_frames))
        phosphates = u.select_atoms("name P31")
        COM = phosphates.center_of_geometry()
        nitrogens_n = u.select_atoms("name N and resname RN SN")
        R_n = nitrogens_n.positions
        phosphates = u.select_atoms("name P31")
        phos = phosphates.positions

        for i in range(0, len(R_n)):
            OUTPUT1.write('{}     '.format(ts.frame))
            OUTPUT1.write('{}     '.format(R_n[i][0]))
            OUTPUT1.write('{}     '.format(R_n[i][1]))
            OUTPUT1.write('{}   \n'.format(R_n[i][2] - COM[2]))

        for i in range(0, len(phos)):
            OUTPUT2.write('{}     '.format(ts.frame))
            OUTPUT2.write('{}     '.format(phos[i][0]))
            OUTPUT2.write('{}     '.format(phos[i][1]))
            OUTPUT2.write('{}   \n'.format(phos[i][2] - COM[2]))

    OUTPUT1.close()
    OUTPUT2.close()
