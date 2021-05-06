#!/bin/sh

for md in WT_EGFR D770_GY  D770_N771insG  N771_P772insN
do
    for run in ../../md/$md/active/*
    do
        echo "0 1" | gmx mindist -f $run/traj_aln.xtc -s $run/prod.tpr -n $run/ke.ndx -od "$md"_"$(basename $run)"_ke.xvg 
    done
    awk '$0 !~ /^[@#]/ {print $2}' "$md"_*_ke.xvg > "$md"_ke.txt
done
