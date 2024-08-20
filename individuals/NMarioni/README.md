# NMarioni PolyAnalysis
The following contains scripts developed, modified, and/or used by Nico Marioni in the Ganesan Polymer Physics group to perform molecular dynamics simulations and post analysis.

**/PolyBuild/** is used to construct .gro and .top gromacs files for polymers of custom length/monomer
 - Gaussian/Multiwfn for charge optimization
 - MKTOP to generate .top for a trimer of desired monomer
 - split_top.py/poly_build.py to disassemble the trimer .top and generate custom polymers
 - README.md inside for more info

**/PolymerCrosslinking/** is used to crosslink a polyamide membrane
 - README.md inside for more info

**/PyAnalysis/** contains useful python post-processing codes
 - **NOTE**: these scripts are provided as is and may require edits to meet your specific needs
 - RDFs/CDFs
 - MSD/Onsager "MSDs" for Correlated Motion/1-D Displacement
   - Base code for Onsager analysis developed by Kara Fong: https://github.com/kdfong/transport-coefficients-MSD
 - Ion Pairing/Cluster Analysis
 - Ion Pair Autocorrelation Functions for extracting Ion Pair Residence Times
 - Water-Oxygen hydrogen bond distributions
 - Static and Dynamic Structure analyses for extracting Polymer Segmental Relaxation Timesscales
 - Pore Size Distributions
 - Extra - Save unwrapped system CoM at every frame for certain analyses, Convert .gro to .xyz for poreblazer

## Usage
/PolyBuild/ and /PolymerCrosslinking/ are explained iwithin the README.md's of those folders

/PyAnalysis/ and Python script general use below

### **Notable** and *Useful* Packages
 - **MDAnalysis** - Used to read in GROMACS .trr/.xtc and .tpr files and do basic analysis
 - **h5py** - Used as a simple and efficient way to save large arrays of data/save memory as a .hdf5 I/O file
 - **networkx** - Used as a simple and efficient way to analyze clusters
 - **multiprocessing/functools** - Used to parallelize python scripts over multiple threads
   - **NOTE - ONLY PARALLELIZES ACROSS 1 NODE**
 - *numpy* - Arrays and array math is mostly handled with numpy
 - *scipy* - Scipy contains several useful functions
 - *memory-profiler* - Useful tool to troubleshoot problems with running out of memory
 - *os* - OS is commonly used to check if files exists and to delete files
   - This is commonly used to separate loading in data (load_TRR() below) from analysis to combat memory problems
 - *matplotlib* - A powerful plotting tool

### Common Functions
 - **load_TRR()** - Load in data from GROMACS .trr/.xtc/.tpr files using MDAnalysis
 - **XXX_analysis()** - Commonly hosts the multiprocessing steps of XXX_calc()
 - **XXX_calc()** - The calculation that is parallelized in XXX_analysis()
    
### Common Arguments
 - **trj_file** - GROMACS .trr/.xtc file name
 - **top_file** - GROMACS .tpr file name
 - **aX_name** - Atom name as shown in .top/.itp
 - **coord** - coordination distance (Angstroms)
 - **t_min** - Start time for analysis (ps)
   - t_min = -1 starts at first frame
 - **t_max** - End time for analysis (ps)
   - t_max = -1 starts at last frame
 - **step** - Frame step-size (step = 1 reads every frame, 2 every other frame, etc)
   - step = -1 is the same as step = 1
 - **nt** - Number of threads used in multiprocessing steps
 - **CoM_check** - If .trr/.xtc/.tpr does not contain all of atoms in the system (commonly done to reduce file size and load times), then the System CoM will be incorrect. CoM_Check = 1 deals with that using CoM.hdf5 genereated from Sys_CoM.py. CoM_Check = 2 can oftin times be used to perform the analysis relative to a different center of mass, e.g., polymer or water.

### Example Code Execution
```
# Computes the RDF and CDF g_NA-CL, n_NA-CL over 20 ns of data from 80-100 ns
python3 ${Insert Path}/analysis_rdf.py md.trr md.tpr NA CL 2.0 1 80000 100000 1 128
    where - script, trj_file, top_file, a1_name, a2_name, rmax, cdf_check, t_min, t_max, step, nt = argv

# Computes the MSD for NA over the entire trajectory, where the system CoM is retrieved from CoM.hdf5
python3 ${Insert Path}/analysis_msd.py unwrap.trr md.tpr NA -1 -1 -1 128 1
    where - script, trj_file, top_file, a_name, t_min, t_max, step, nt, CoM_check = argv
```

## **Necessary** and *Useful* Software
**Gromacs** - Version 2020.5

**Gaussian** - g16 is available through TACC

**Multiwfn** - Version 3.7, Binary, Linux, No GUI
 - Link provided within /PolyBuild/

**MKTOP** - (Modified) Version 2.2.1
 - Github repository and necessary modifications are listed within /PolyBuild/

**Python** - Version 3.8.10
```
# Current Python Packages: python3 -m pip list
Package         Version
--------------- -------
biopython           1.83
contourpy           1.1.1
cycler              0.12.1
fasteners           0.19
fonttools           4.53.1
GridDataFormats     1.0.1
gsd                 3.2.1
h5py                3.11.0
importlib-resources 6.4.0
intermol            0.1.0.dev0
joblib              1.4.2
kiwisolver          1.4.5
matplotlib          3.7.5
MDAnalysis          2.4.3
memory-profiler     0.61.0
mmtf-python         1.1.3
mrcfile             1.5.2
msgpack             1.0.8
networkx            3.1
numpy               1.24.4
packaging           24.1
ParmEd              4.2.2
pillow              10.4.0
pip                 21.1.1
psutil              6.0.0
pyparsing           3.1.2
python-dateutil     2.9.0.post0
scipy               1.10.1
setuptools          56.0.0
six                 1.16.0
threadpoolctl       3.5.0
tqdm                4.66.4
zipp                3.19.2
```

*PoreBlazer* - Version 4.0
 - https://github.com/SarkisovGitHub/PoreBlazer
 - Great tool for calculating pore size distributions
 - /PyAnalysis/analysis_PSD.py was developed using similar methods

*InterMol* - Version 0.1.2
 - https://github.com/shirtsgroup/InterMol
 - Great tool for converting between different simulation software
 - I have used this frequently to convert from GROMACS to LAMMPS