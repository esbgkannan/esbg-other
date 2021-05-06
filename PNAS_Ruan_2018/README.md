MD Data Files for Manuscript "Altered conformational landscape and dimerization dependency underpins the activation of EGFR by αC–β4 loop insertion mutations".
===============================================================================================================================================================

The repository contains
-----------------------

* MD trajectories for WT EGFR, D770_N771insG, N771_P772insN, D770>GY, L704N, and L760R.
* Umbrella sampling data points for WT EGFR, D770_N771insG, N771_P772insN, and D770>GY.
* shell/python code to generate key figures in the paper.

Package dependencies
--------------------
Please ensure the following software is installed:

- `Gromacs 5.0+ <http://www.gromacs.org>`
- `NumPy <http://www.numpy.org>`
- `pyEMMA <http://emma-project.org/latest/>`
- `matplotlib <https://matplotlib.org>`


MD trajectories
---------------

The full MD trajectories are around 37GB and are accessible through [this dropbox link](https://www.dropbox.com/sh/zw8w6wnptuu2d7g/AADekNzCm8XZwc-s5LoytKV8a?dl=0).
Note the file hierarchy should be stored under `md/` folder for the shell/python script files to work.
Each MD file folder contains the following files:

* traj_aln.xtc -- Gromacs trajectory file with non-protein atoms removed and protein aligned
* prod.tpr -- Gromacs topology file to run the simulation
* confout.gro -- Structure file of the protein in gro format.
* index.ndx -- Gromacs index file for the simulation. Note that non-protein atoms are removed in traj_aln.xtc.


Umbrella sampling data
----------------------
The umbrella sampling data (US) is located in `us/`. CV files for three independent replicates are included for each mutant.
The column of each cv files corresponds to `shapshot number`, `CV1` and `CV2`.

Code to generate key figures in the paper
-----------------------------------------
The `code/` directory contains shell/python code to generate fig 5a, fig 6a and fig 7a in the paper

#### Fig 5a
Make sure you have downloaded the MD trajectories from [the dropbox link](https://www.dropbox.com/sh/zw8w6wnptuu2d7g/AADekNzCm8XZwc-s5LoytKV8a?dl=0)
and store them in the `md/` folder with exactly the same directory tree.

    1. Generate the RMSD for C-helix of EGFR MD simulation. In a shell environment, run

        # cd code/fig5a/chelix_rmsd/
        # sh get_chelix_rmsd.sh

    2. Generate K745-E762 distance of EGFR MD simulation. In a shell environment, run
        
        # cd ../ke_distance/
        # sh get_ke_distance.sh

    3. Generate the 2D distribution plot. Make sure you have numpy and matplotlib package installed. In a shell environment, run
        
        # cd ../make_plot/
        # python plot.py 


#### Fig 6a
Make sure you have downloaded the MD trajectories from [the dropbox link](https://www.dropbox.com/sh/zw8w6wnptuu2d7g/AADekNzCm8XZwc-s5LoytKV8a?dl=0)
and store them in the `md/` folder with exactly the same directory tree.

    1. Generate the R776-A767 capping distance. In a shell environment, run

        # cd code/fig6a/r776-a767_distance/
        # sh get_capping_distance.sh

    3. Generate the 2D distribution plot. Make sure you have numpy and matplotlib package installed. In a shell environment, run
        
        # cd ../make_plot/
        # python plot.py 

#### Fig 7a
The calculation of free energy landscape (FEL) requires the [pyEMMA package](http://emma-project.org/latest/).
The python script `free_energy.py` reads all the US data and performs weighted histogram analysis method (WHAM) calculation.
Follow this procedure to generate an FEL for WT EGFR.

    1. Generate a US data file:

        # cd code/fig7a/
        # ls -v ../../us/WT/us_data1/cv_*txt > us_data1_cv.txt

    2. Perform WHAM calculation and plot the FEL. If the US data file has a differnet name other than us_data1_cv.txt, modify line 10 of draw_pmf.py accordingly.
        
        # python free_energy.py
        # python draw_pmf.py


