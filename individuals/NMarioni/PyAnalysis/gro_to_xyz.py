#!/usr/bin/python

#Converts .gro file to .xyz for use in poreblazer

from sys import argv
scipt, filename = argv
# filename = name of .gro file (not including '.gro' to be converted to {$filename}.xyz)

# check defines the residue names that are removed from the .gro file
check = ['LI', 'NA', 'K', 'CL', 'SOL']; brk = 0

with open(filename+'.gro','r') as f:
    gro = f.readlines()

write_file = open(filename+'.xyz', 'w')

write_file.write(gro[1])
write_file.write('XYZ\n')

count = 0
for line in gro[2:len(gro)-1]:
    line = line.split()

    for c in check:
        if line[0][len(line[0])-len(c):] == c:
            brk = 1; break
    if brk == 1:
        break

    if len(line) == 9:
        atom = line[1][0]
        x = float(line[3])*10; y = float(line[4])*10; z = float(line[5])*10
    elif len(line) == 8:
        atom = line[1][0]
        x = float(line[2])*10; y = float(line[3])*10; z = float(line[4])*10

    write_file.write('{:s} {:f} {:f} {:f}\n'.format(atom, x, y, z))
    count += 1

write_file.write('\n')

write_file.close()
box = str(round(float(gro[len(gro)-1].split()[0])*10, 4))

with open(filename+'.xyz','r') as f:
    gro = f.readlines()

write_file = open(filename+'.xyz', 'w')
write_file.write(str(count)+'\n')
for line in gro[1:]:
    write_file.write(line)
write_file.close()

write_file = open('input.dat', 'w')
write_file.write(filename+'.xyz\n')
write_file.write('{:15s}{:15s}{:15s}\n'.format(box,box,box))
write_file.write('90             90             90\n\n')
write_file.close()
