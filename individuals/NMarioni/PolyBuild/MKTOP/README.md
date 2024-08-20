# MKTOP
**create_top** - shell script to 1) convert monomer head/mid/tail .pdb's to .gro's and 2) use MKTOP to create monomer trimer .top for use with GROMACS

**mktop.pl** - slightly modified MKTOP script to generate a .top topology file from a given pdb
 - **NOTE**: Code is specifically designed to work with .top files created by an edited version of MKTOP. All edits are minor changes to the formatting of the generated .top file. MKTOP can be found here: https://github.com/aar2163/MKTOP Edits are outlined below.
 - **NOTE**: Line 3 $gromacs_dir="{PATH}" should path to your gromacs directory
 - **NOTE**: Lines 640, $ffdir = "amber03.ff", and 644, $ffdir = "oplsaa.ff", should specify your forcefield folder of interst
 - **NOTE**: Lines 2858, 2870, 2883, 2907, 2922 should have only 1 space between 1 and improper_X_X_X_X, e.g.:
   <pre> $improper[$I]="$k[0] $k[1] $i $j 1 improper_O_C_X_Y" **NOT** $improper[$I]="$k[0] $k[1] $i $j 1    improper_O_C_X_Y" </pre>
 - **NOTE**: Lines 3071-3118 should be replaced by the following code:

   ```
   for ($i = 1; $i <= $linha; $i++) {
    $Z=$Z[$i];
       printf OUT "   %3d", $i;
       print OUT "   $type[$i]      1     $res[$i]";
       $atom=$mol[$Z];
       printf OUT "  %4s", $atom;
       printf OUT "     %2d", $group[$i];
       printf OUT "   % 16.15f", $q[$i];
       printf OUT "    %7.4f", $mass[$Z];
       print OUT "\n";
          $count[$Z]++;
      }
   
   
   print "\nWriting all the stuff...\n";
   
   print OUT "\n\n";
   
    print OUT "\n\n[ bonds ]\n";
    for ($i=1; $i < $B; $i++) {
     @list = split (/ /, $BOND[$i]);
     printf OUT"%4d %4d %d", $list[0], $list[1], $list[2];
     #printf OUT "%6.3f  ", ($list[3]/10);
     #print OUT "$list[4]";
     print OUT "\n";
    }
   
    print OUT "\n\n[ angles ]\n";
    for ($i=1; $i < $A; $i++) {
     @list = split (/ /, $angle[$i]);
     printf OUT "%4d %4d %4d %d", $list[0], $list[1], $list[2], $list[3];
     #printf OUT "%6.3f  ", ($list[4]);
     #print OUT "$list[5]";
     print OUT "\n";
    }
   
    print OUT "\n\n[ dihedrals ]\n";
    for ($i=1; $i < $D; $i++) {
     @list = split(/ /, $dihedral[$i]);
    # if ($list[5] != 0 or $list[6] != 0 or $list[7] != 0 or $list[8] != 0) {
      printf OUT "%4d %4d %4d %4d %d", $list[0], $list[1], $list[2], $list[3], $list[4];
      print OUT "\n";
    # }
    }
    print OUT "\n\n[ dihedrals ]\n";
    for ($i=1; $i < $I; $i++) {
     @list = split(/ /, $improper[$i]);
      printf OUT "%4d %4d %4d %4d %d %s", $list[0], $list[1], $list[2], $list[3], $list[4], $list[5];
      print OUT "\n";
    }
    print OUT "\n\n[ pairs ]\n";
    for ($i=1; $i < $P; $i++) {
     @list = split(" ", $pair[$i]);
      printf OUT"%4d %4d %d", $list[0], $list[1], $list[2];
      print OUT "\n";
    }
   ```

**split_top.py** - python script that disassembles the trimer .top and monomer head/mid/tail .gro files for use in poly_build
 - **NOTE**: It is important to use pdb2gmx to create .gro files and the above modified MKTOP to create .top. Any changes in formatting could lead to errors
 - **NOTE**: There are commented out print statements at the end of many blocks of code. This can be useful for validation/bug fixing

## Example
/DVB and /SPS are examples for Divinylbenzene and Sulfonated Styrene
 - **NOTE**: .pdb files should be created in Avogadro or similar software and inserted into folder
   - ***IMPORTANT***: Be careful when creating your trimer/head/mid/tail .pdbs. Numbering systems should be consistent between units with head/tail caps numbered at the end.
 - **NOTE**: _charge.dat file is created using knwon charges (i.e., OPLS charges), or obtained from charge optimization/RESP charges using Gaussian/Multiwfn, etc
 - **NOTE**: .gro and .top files are created using the create_top shell script
 - **NOTE**: top.hdf5 is created by split_top.py for use by poly_build






**MKTOP README BELOW - See the github provided above to access the code**

*** PLEASE CITE THIS WORK ***

Ribeiro, A.A.S.T.; Horta, B.A.C.; de Alencastro, R.B.  J. Braz. Chem. Soc., Vol. 19, No. 7, 1433-1435, 2008.

*************************************

MKTOP version 2.2.1

MKTOP is a Perl Script designed to build topologies for GROMACS (4.5.x)

If you are using Linux/Unix, the command line would be:

mktop_2.2.1.pl -i input.pdb -c charges.txt -o topology.top -ff opls

*** INPUT FORMAT *****

From time to time someone contacts me with problems running MKTOP. While I
definitively appreciate any feedback regarding possible bugs, most reported
errors are actually improper input formating. From version 2.0 on, element
symbol on columns 77-78 has been obligatory. IT SHOULD BE RIGHT-JUSTIFIED AND
UPPER-CASE.

***CHARGES ***
The (optional) input file charges.txt contains atomic charges that you must calculate by some method (we recommend using RESP).
This file must have the following syntax:
(number of atom) (charge)

So, if you decide to perform a simulation of sodium chloride in vacuum, the file would be:

1 1.0000
2 -1.000

If you do not provide a charges.txt file we will assume zero charges for all
atoms!!!

The numbering in the files input.pdb and charges.txt must be the same!!!

***FORCEFIELD ***

MKTOP now supports both AMBER03 and OPLS-AA forcefields, so you need to set
the -ff option either to amber or opls

***CONNECTIVITY****
If your PDB file has connectivity information (lines starting with CONECT) you
can use the -conect option set to yes. MKTOP will use this info to determine
bonds, angles and dihedrals. If you do not have this info, you have to set
-conect equal to no so the script will get this info from the distances in
your PDB. 

In the former case pay extra attention to your output. MKTOP will
use distances between atoms to determine bonded pairs. If some of the
distances in your PDB fall outside of the "bonded range", EVERYTHING WILL
FAIL.

***TIPS ***
Do not forget to change the gromacs_dir line in the beginning of the script to your actual path.

If you have never executed anything on your command line, you are supposed to change the permissions of the file mktop_2.0.pl to make it executable:

chmod u+rxw mktop.pl

If you have any questions, send them to andre@aribeiro.net.br

