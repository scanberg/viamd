import math

file = open("test_4_points.pdb", "w")

file.write("CRYST1   10.000   10.000   10.000  90.00  90.00  90.00    \n")

frames = 2000
dt = 0.4
for i in range(0, frames):
    file.write("MODEL     %i\n" % (i))
    for j in range(4):
        #offset = j * 0.35
        x = 2.5 + (j < 2) * 5.0
        y = 2.5 + (j % 2 == 0) * 5.0
        z = 5.0 + math.sin(dt * i * (1.0 + j * 0.03)) * 5.0
        occ = 1.0
        tmp = 0.0
        file.write("ATOM      %i  C   RES     1    %8.3f%8.3f%8.3f  1.00  0.00           C  \n" % (j + 1, x, y, z))
    file.write("ENDMDL\n")