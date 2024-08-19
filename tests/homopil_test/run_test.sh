#! /bin/sh
#
# run_test.sh
# Copyright (C) 2019 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


SRC="../../src"

python $SRC/utanalysis.py --top topol.tpr --trj traj.xtc --end 20 --func homopoly --nt 2 && \
python $SRC/utanalysis.py --top topol.tpr --trj traj.xtc --end 20 --func aacf --nt 2 && \
python $SRC/utanalysis.py --top topol.tpr --trj traj.xtc --end 20 --func hopping --nt 2 || exit 1

mkdir -p test_dat
mv *.dat test_dat

diff -Nur test_dat reference
DIFF_RESULT=$?

if [ $DIFF_RESULT = 0 ]; then
    rm -rf test_dat
    rm -rf utanalysis.log
fi

exit $DIFF_RESULT
