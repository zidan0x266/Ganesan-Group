import MDAnalysis as mda
import numpy as np
uta = mda.Universe('../cg_topol.tpr', '../cg_traj.xtc')
aP = uta.select_atoms("type PF")
for ts in uta.trajectory:
    if ts.frame % 50 == 0:
        print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, uta.trajectory.time))
        cell = uta.trajectory.ts.dimensions
        ifm1 = [8.0]; ifn1 = [8.0, 26.9]
        bkm1 = [26.9, 53.1]; bkn1 = [76.9, 100.1]
        ifm2 = [53.1, 76.9]; ifn2 = [100.1]
        # local range processing
        Ainif1 = []; Ainif2 = []; Ainif3 = []; Ainif4 = []  # Atoms in interface
        Ainbk1 = []; Ainbk2 = []  # Atoms in bulk
        for i in range(len(aP)):
            # interface
            if (aP.positions[:,2][i] < ifm1[0]):
                Ainif1.append(i)
            if (aP.positions[:,2][i] >= ifm2[0] and aP.positions[:,2][i] < ifm2[1]):
                Ainif2.append(i)
            if (aP.positions[:,2][i] >= ifn1[0] and aP.positions[:,2][i] < ifn1[1]):
                Ainif3.append(i)
            if (aP.positions[:,2][i] >= ifn2[0]):
                Ainif4.append(i)
            # bulk
            if (aP.positions[:,2][i] >= bkm1[0] and aP.positions[:,2][i] < bkm1[1]):
                Ainbk1.append(i)
            if (aP.positions[:,2][i] >= bkn1[0] and aP.positions[:,2][i] < bkn1[1]):
                Ainbk2.append(i)
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
    if ts.frame == 1000:
        break
print("Job finished!")
