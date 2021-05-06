#!/bin/sh

for md in WT_EGFR D770_GY  D770_N771insG  N771_P772insN
do
    for run in ../../md/$md/active/*
    do
        echo "0 0" | gmx rms -f $run/traj_aln.xtc -s $run/prod.tpr -n $run/chelix.ndx -o "$md"_"$(basename $run)"_chelix_rmsd.xvg 
    done
    awk '$0 !~ /^[@#]/ {print $2}' "$md"_*_chelix_rmsd.xvg > "$md"_chelix_rmsd.txt
done
