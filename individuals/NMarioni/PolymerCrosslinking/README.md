# PolymerCrosslinking
**cross_link.py** is a python script developed to take a GROMACS .gro and .top file and perform a crosslink step based on the current configuration
 - **NOTE**: This code is provided as is, and is not designed for general use

**/mdp/** contains useful supporting files, including the files necessary to crosslink a TMC + MPD polyamide membrane
 - **NOTE**: /mdp/ contains the builder_io.py script, developed by Jakud Krajniak and Zidan Zhang, modified by Nico Marioni, used in crosslink.py and mdp/sys_build.py
 - **NOTE**: mdp/sys_build.py is used to create mdp/polyamide.top

 **run.sh** is a bash script designed to crosslink a polyamide membrane using GROMACS and cross_link.py
  - see comments in both scripts for more information