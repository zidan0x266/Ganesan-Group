# Gaussian/Multiwfn
**fchk.sh** - converts .chk files from Gaussian (https://gaussian.com/) into .fchk files for use by Multiwfn

**gXX_test** - shell script to run Gaussian

**Multiwfn** - executable file used to obtain RESP charges from Gaussian optimized structure/charges
 - **NOTE**: Multiwfn can be found here: http://sobereva.com/multiwfn/

**settings.ini** - Multiwfn settings file
 - **NOTE**: Line 92 gaupath, Line 93 cubegenpath must path to the gXX and cubegen file in your Gaussian install
 - **NOTE**: Line 95 formchkpath can optionally be pathed to formchk in your Gaussian install to allow Multiwfn to convert .chk to .fchk

## Example
/DVB/ and /SPS/ are examples for Divinylbenzene and Sulfonated Styrene (CR61)
 - **NOTE**: *.com files are supplied to run Gaussian (see below), .log/.chk are created by Gaussian, .fchk is made with formchk for use by Multiwfn
 - **NOTE**: *_chg.txt/_con.txt are supplied based on your needs for conformations/charge constraints for use with Multiwfn (see below)
 - **NOTE**: *_charges is the RESP charges output by Multiwfn. This becomes the *_charge.dat files in /MKTOP/

## Gaussian: Strucutre/Charge Optimization
Avogadro (https://avogadro.cc/) gives gaussion *.com files styled as shown below (>Extensions>Gaussian for template)

```
%NProcShared=128                                  #Number of threads
%mem=200GB                                        #Memory
%Chk=DVB_iso.chk                                  #Name of .chk file to be made
#n B3LYP/6-311+g** Opt                            #Type of optimization

 DVB_iso                                          #File Name

0 1                                               #First = Total Charge, Second = Multiplicity
C         -3.51176       -0.14610        2.98652  #Starting file coords in xyz
```

## Multiwfn: RESP Charges
It can be useful to compute charges averaged over multipled conformations (ex: isotactic and syndiotactic conformations). If used, create *_con.txt (example below)

```
LGB_syn.fchk 0.5    # First conformation .fchk file with weight (0.5 means it will make up 50% of the final charges)
LGB_iso.fchk 0.5    # Second conformation .fchk with weight (weights should add to 1.0)
```

You may need to apply charge constraints (ex: for a trimer you want each unit (head/mid/tail) to be constrained to their respective charge). If used, create *_chg.txt (example below)

```
1-26 0    #First monomer (head unit) with charge (0, -1, +1, etc)
27-51 0   #Second monomer (mid unit) with charge (third monomer (tail unit) is assumed, i.e. all non-specified atom numbers will be constrained to the remaining charge to achieve total charge from .com file above))
```

### Running Multiwfn
 - **NOTE**: IDK of any way to automate this
 - **NOTE**: Copy the charges from the printed screen (as seen in *_charges file in example folder)

```
../Multiwfn ${file_name}_syn.fchk
7   # Population analysis and atomic charges
18  # Restrained ElectroStatic Potential (RESP) atomic charge

# If using Multiple Conformations
-1  # Load list of conformer and weights from external file
${file_name}_con.txt

# If using Charge Constraints
6   # Set charge constraint in fitting, current: No constraint
1   # Load charge constraint setting from external plain text file
${file_name}_chg.txt

2 # Start one-stage ESP fitting calculation with constraints
q
```
