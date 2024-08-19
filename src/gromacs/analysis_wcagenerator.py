#!/usr/bin/env python
#
#  Python code to calculate VdW tabulate potential for Gromacs
#
#  Copyright 2018-2021, Zidan Zhang <zhangzidan@gmail.com>
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  wcagenerator.py
#
#     Usage:
#     python wcaGenerator.py GROUPA GROUPB SIGB EPSB FUN
#     python wcaGenerator.py A B 1.0 1.0 2.5 1


import math
import string
import argparse

parser = argparse.ArgumentParser('Generator for WCA potentials')
parser.add_argument('groupa', action='store',
                    help='Name of group 1')
parser.add_argument('groupb', action='store',
                    help='Name of group 2')
parser.add_argument('sigma', action='store',
                    type=float,
                    help='LJ Sigma')
parser.add_argument('epsilon', action='store',
                    type=float,
                    help='LJ Epsilon')
parser.add_argument('rcut', action='store',
                    type=float,
                    help='cutoff for shift')
parser.add_argument('fun', action='store',
                    type=int, default=1,
                    help='Function type for LJ interactions')
# Parse arguments
args = parser.parse_args()
# Call function

def ljwcat1(groupa, groupb, sigma, epsilon, rcut):

    # Table Parameters
    xmax = 3.5               # in nm
    dx = 0.002               # in nm
    bmax = int(xmax/dx)

    # Create the table
    fname="table_" + groupa.strip() + "_" + groupb.strip() + ".xvg"
    f = open(fname,'w')

    # Set output format
    fmt = "%14.5e" * 7 + "\n"

    # Generate the table
    for i in range(bmax+1):
        x = i * dx
        if (x == 0) :
            uwca = 4*epsilon*((sigma/(0.5*dx))**12-(sigma/(0.5*dx))**6) + epsilon
            f.write( fmt % (0.0,0.0,0.0,0.0,0.0,uwca,0.0) )
        else:
            if (x <= rcut**(1.0/6.0)*sigma):
                uwca = 4*epsilon*((sigma/x)**12-(sigma/x)**6) + epsilon
                if (i == 1):
                    uwcam1 = 4*epsilon*((sigma/(0.5*dx))**12-(sigma/(0.5*dx))**6) + epsilon
                else:
                    xm1 = x - dx
                    uwcam1 = 4*epsilon*((sigma/xm1)**12-(sigma/xm1)**6) + epsilon
                xp1 = x + dx
                uwcap1 = 4*epsilon*((sigma/xp1)**12-(sigma/xp1)**6) + epsilon
                fwca = (uwcam1 - uwcap1)/2*dx
                f.write( fmt %
                         (x,0.0,0.0,0.0,0.0,uwca,fwca) )
            else:
                f.write( fmt %
                         (x,0.0,0.0,0.0,0.0,0.0,0.0) )

    f.close()

def ljwcat2(groupa, groupb, sigma, epsilon, rcut):

    # Table Parameters
    xmax = 3.5               # in nm
    dx = 0.002               # in nm
    bmax = int(xmax/dx)

    rc = rcut**(1.0/6.0)*sigma

    # Create the table
    fname="table_" + groupa.strip() + "_" + groupb.strip() + ".xvg"
    f = open(fname,'w')

    # Set output format
    fmt = "%14.5e" * 7 + "\n"

    # Generate the table
    for i in range(bmax+1):
        x = i * dx
        if (x == 0) :
            uwca = 4*epsilon*((sigma/(0.9*dx))**12-(sigma/0.9*dx)**6-(sigma/rc)**12+(sigma/rc)**6)
            f.write( fmt % (0.0,0.0,0.0,0.0,0.0,uwca,0.0) )
        else:
            if (x <= rc):
                uwca = 4*epsilon*((sigma/x)**12-(sigma/x)**6-(sigma/rc)**12+(sigma/rc)**6)
                if (i == 1):
                    uwcam1 = 4*epsilon*((sigma/(0.9*dx))**12-(sigma/0.9*dx)**6-(sigma/rc)**12+(sigma/rc)**6)
                else:
                    xm1 = x - dx
                    uwcam1 = 4*epsilon*((sigma/xm1)**12-(sigma/xm1)**6-(sigma/rc)**12+(sigma/rc)**6)
                xp1 = x + dx
                uwcap1 = 4*epsilon*((sigma/xp1)**12-(sigma/xp1)**6-(sigma/rc)**12+(sigma/rc)**6)
                fwca = (uwcam1 - uwcap1)/2*dx
                f.write( fmt %
                         (x,0.0,0.0,0.0,0.0,uwca,fwca) )
            else:
                f.write( fmt %
                         (x,0.0,0.0,0.0,0.0,0.0,0.0) )

    f.close()

def purelj(groupa, groupb, sigma, epsilon, rcut):

    # Table Parameters
    xmax = 3.5               # in nm
    dx = 0.002               # in nm
    bmax = int(xmax/dx)

    # Create the table
    fname="table_" + groupa.strip() + "_" + groupb.strip() + ".xvg"
    f = open(fname,'w')

    # Set output format
    fmt = "%14.5e" * 7 + "\n"

    # Generate the table
    for i in range(bmax+1):
        x = i * dx
        if (x == 0) :
            uwca = 4*epsilon*((sigma/(0.5*dx))**12-(sigma/(0.5*dx))**6)
            f.write( fmt % (0.0,0.0,0.0,0.0,0.0,uwca,0.0) )
        else:
            uwca = 4*epsilon*((sigma/x)**12-(sigma/x)**6)
            if (i == 1):
                uwcam1 = 4*epsilon*((sigma/(0.5*dx))**12-(sigma/(0.5*dx))**6)
            else:
                xm1 = x - dx
                uwcam1 = 4*epsilon*((sigma/xm1)**12-(sigma/xm1)**6)
            xp1 = x + dx
            uwcap1 = 4*epsilon*((sigma/xp1)**12-(sigma/xp1)**6)
            fwca = (uwcam1 - uwcap1)/2*dx
            f.write( fmt %
                     (x,0.0,0.0,0.0,0.0,uwca,fwca) )

    f.close()


def main():
    if (args.fun == 1):
        ljwcat1(args.groupa, args.groupb, args.sigma,
                args.epsilon, args.rcut)
    elif (args.fun == 2):
        ljwcat2(args.groupa, args.groupb, args.sigma,
                args.epsilon, args.rcut)
    elif (args.fun == 3):
        purelj(args.groupa, args.groupb, args.sigma,
                args.epsilon, args.rcut)


if __name__ == '__main__':
    main()