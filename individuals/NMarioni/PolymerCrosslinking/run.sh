#!/bin/bash

# This script will create a system with 690 MPD and 460 TMC monomers, and then crosslink it to 90%

nMPD=690
tMPD=$(($nMPD-1))
nTMC=460

threads=48

# polyamide.top can be created within the mdp folder using sys_build.py
cp mdp/polyamide.top ./topol.top
seed=$(( RANDOM % 10000 ))
echo $seed
gmx insert-molecules -f mdp/MPD.gro -ci mdp/MPD.gro -nmol $tMPD -try 1000 -box 5.04 5.04 10.08 -seed $seed -o initial.gro
gmx insert-molecules -f initial.gro -ci mdp/TMC.gro -nmol $nTMC -try 1000 -seed $seed -o solv.gro -scale 0.40
rm initial.gro



mkdir -p equil

# Short energy minimzation
gmx grompp -f mdp/em.mdp -c solv.gro -p topol.top -o em
ibrun mdrun_mpi -deffnm em
rm *# step*
mkdir -p equil/em
mv em* equil/em
rm equil/em/*.trr equil/em/*.xtc

# The system will be crosslinked with non-pbc in the z-direction, but pbc in the x and y directions
# Vacuum is added to ensure crosslinks cannot be accidentally formed across the z boundary
printf '0\n' | gmx trjconv -f equil/em/em.gro -s equil/em/em.tpr -o save.gro -pbc res
gmx editconf -f save.gro -o save.gro -box 5.04 5.04 10.5
cp topol.top save.top



# Crosslink the membrane
cp save.top topol.top
cp save.gro eq.gro
rm xlink_tracker.txt
i=1
j=0
while :
do
    echo "Cycle $i"
    # Perform a crosslink step
    python3 cross_link.py eq.gro eq.gro topol.top topol.top
    rm *_

    # If no new crosslinks are formed for 50 consequetive NVT simulations, or the Crosslinked Density reached 90% -> Crosslinking is complete -> end process
    if [ -f "xlink_tracker.txt" ]; then
        tracker="$(cat xlink_tracker.txt)"
        echo "Tracker $tracker"
        if [ "$tracker" -gt "49" ]; then
            echo "Done Crosslinking"
            break
        fi
    else
        echo "Still Crosslinking"
    fi

    # Short energy minimzation
    gmx grompp -f mdp/xlink_em.mdp -c eq.gro -p topol.top -o eq
    ibrun mdrun_mpi -deffnm eq
    rm *# step*
    
    # 10 ps NVT simulation before the next crosslink step
    gmx grompp -f mdp/xlink_NVT.mdp -c eq.gro -p topol.top -o eq
    ibrun -np 1 mdrun_mpi -deffnm eq
    rm *#

    # If NVT simulation fails -> start over
    if ! grep -q 'Finished mdrun' eq.log; then
        i=0
        cp save.top topol.top
        cp save.gro eq.gro
        ((j++))
        echo "System Restarted $j times"
    fi

    # If start over too many times -> end script -> System Failed
    if [ "$j" -gt "4" ]; then
        echo "System Failed"
        exit
    fi

    ((i++))
done
mkdir -p equil/xlink
mv eq* equil/xlink
rm equil/xlink/*.trr equil/xlink/*.xtc
