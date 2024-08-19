#! /bin/sh
#
# run_tests.sh
# Copyright (C) 2019 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

oldpwd="`pwd`"

for test_dir in *_test; do
    cd $test_dir
    echo "running ${test_dir}"
    ./run_test.sh 
    RES=$?
    if [ $RES != 0 ]; then
        exit $RES
    fi
    cd $oldpwd
done
