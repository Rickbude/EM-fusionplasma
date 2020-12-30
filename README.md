## Introduction

This directory contains the MATLAB code that was used to generate the 
results that can be found in the paper:

**'Accelerating Simulations of Electromagnetic Waves in Hot, Magnetized Fusion Plasmas'**

by R.H.S. Budé, D. Van Eester, J. van Dijk, R.J.E. Jaspers and A.B. Smolders,
accepted for publication in PPCF (Plasma Physics and Controlled Fusion). 
doi: https://doi.org/10.1088/1361-6587/abd619

## License

The code has been made available under the GNU GENERAL PUBLIC LICENSE version 3,
see the file LICENSE.md for details.

## Using the code

- Create an input file in the subdirectory input_files. You can use one of the
  existing files as an example. The results in the paper have been generated
  with the input file "JET_case_paper.m".
- "Activate" the input file, by running load_input.m. This script will ask for
  the name of one of the input files in the command window. The name of the
  input file will be stored in input_file.mat, and it will subsequently get
  loaded on every run of the code.
- To generate the same plots as in the paper, run the scripts in the
  subdirectory "paper_plots". The electric field calculation happens in
  `hot_plasma_functions/hp_calc_calc_E_truncated_taylor.m`,
  `hot_plasma_functions/hp_calc_E_truncated_fit.m` and
  `hot_plasma_functions/hp_calc_E_all_orders.m`. See the files in "paper_plots",
  and `hot_plasma_not_uniform_1_5_D.m` for examples of how to set up the
  simulations, and of what can be done with the results.  

## Comments, Bugs

If you have any questions or comments about this code, please contact the author
via r.h.s.bude@tue.nl (Rick Budé), or contact Jan van Dijk (j.v.dijk@tue.nl),
Roger Jaspers (r.j.e.jaspers@tue.nl), Bart Smolders (a.b.smolders@tue.nl) or
Dirk Van Eester (d.van.eester@fz-juelich.de)
