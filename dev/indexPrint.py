import MDAnalysis as mda
import numpy as np
uta = mda.Universe('cg_topol.tpr', 'cg_traj.xtc')
aP = uta.select_atoms("type PF")
for ts in uta.trajectory:
    if ts.frame % 100 == 0:
        print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, uta.trajectory.time))
        cell = uta.trajectory.ts.dimensions
        if1 = [23.3, 31.3]; if2 = [53.0, 59.8]; if3 = [75.2, 83.0]; if4 = [100.9]
        bk1 = [31.3, 53.0]; bk2 = [83.0, 100.9]
        # local range processing
        Ainif1 = []; Ainif2 = []; Ainif3 = []; Ainif4 = []  # Atoms in interface
        Ainbk1 = []; Ainbk2 = []  # Atoms in bulk
        for i in range(len(aP)):
            # interface
            if (aP.positions[:,2][i] >= if1[0] and aP.positions[:,2][i] < if1[1]):
                Ainif1.append(aP.indices[i] + 1)
            if (aP.positions[:,2][i] >= if2[0] and aP.positions[:,2][i] < if2[1]):
                Ainif2.append(aP.indices[i] + 1)
            if (aP.positions[:,2][i] >= if3[0] and aP.positions[:,2][i] < if3[1]):
                Ainif3.append(aP.indices[i] + 1)
            if (aP.positions[:,2][i] >= if4[0]):
                Ainif4.append(aP.indices[i] + 1)
            # bulk
            if (aP.positions[:,2][i] >= bk1[0] and aP.positions[:,2][i] < bk1[1]):
                Ainbk1.append(aP.indices[i] + 1)
            if (aP.positions[:,2][i] >= bk2[0] and aP.positions[:,2][i] < bk2[1]):
                Ainbk2.append(aP.indices[i] + 1)
        Ainif = sorted(Ainif1 + Ainif2 + Ainif3 + Ainif4)
        Ainbk = sorted(Ainbk1 + Ainbk2)
        Ainto = sorted(Ainif + Ainbk)
        with open('interface.ndx', 'a+') as anaout:
            print("[ Frame: {0:5d}, Time: {1:8.3f} ps ]".format(ts.frame, uta.trajectory.time), file=anaout)
            for i in range(len(Ainif)):
                if (i + 1) % 8 == 0:
                    print('{:4d}'.format(Ainif[i]), file=anaout, end = '\n')
                elif i == len(Ainif) - 1:
                    print('{:4d}'.format(Ainif[i]), file=anaout, end = '\n')
                else:
                    print('{:4d}'.format(Ainif[i]), file=anaout, end = ' ')
        with open('bulk.ndx', 'a+') as anaout:
            print("[ Frame: {0:5d}, Time: {1:8.3f} ps ]".format(ts.frame, uta.trajectory.time), file=anaout)
            for i in range(len(Ainbk)):
                if (i + 1) % 8 == 0:
                    print('{:4d}'.format(Ainbk[i]), file=anaout, end = '\n')
                elif i == len(Ainbk) - 1:
                    print('{:4d}'.format(Ainbk[i]), file=anaout, end = '\n')
                else:
                    print('{:4d}'.format(Ainbk[i]), file=anaout, end = ' ')
        with open('total.ndx', 'a+') as anaout:
            print("[ Frame: {0:5d}, Time: {1:8.3f} ps ]".format(ts.frame, uta.trajectory.time), file=anaout)
            for i in range(len(Ainto)):
                if (i + 1) % 8 == 0:
                    print('{:4d}'.format(Ainto[i]), file=anaout, end = '\n')
                elif i == len(Ainto) - 1:
                    print('{:4d}'.format(Ainto[i]), file=anaout, end = '\n')
                else:
                    print('{:4d}'.format(Ainto[i]), file=anaout, end = ' ')
    if ts.frame == 10000:
        print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, uta.trajectory.time))
    #    break
print("Job finished!")
