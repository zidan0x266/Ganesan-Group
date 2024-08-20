#!/bin/bash
poly_names='DVB SPS' # list of names of polymer units to be used, separated by spaces
HT_name='DVB'        # default 'None'; name of polymer unit to be used as the Head/Tail unit of the Copolymer
num_units='5 5'    # list of numbers of polymer units corresponding to each name in poly_names; e.g. poly_name = 'PEG LGA' num_units = '10 5' means a 15mer of 10 PEG, 5 LGA
charge_scale=0.85    # charge scale to be applied; e.g. charge_scale = 0.85 => 85% charge scale
poly_type=0          # type of polymer; 0 -> Random Statisical Polymer, 1 -> Block Copolymer (will ignore HT_name)

Conc_list='0.04 0.10 0.20 0.50 1.00'
Sample_list='1'

Conc_list='1.00'
#Sample_list='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40'

for i in $Conc_list; do
    drc=conc_$i
    #drc=$i

    for j in $Sample_list; do
        python3 poly_build.py "${poly_names}" "${HT_name}" "${num_units}" ${charge_scale} $drc $j $poly_type
    done
done
