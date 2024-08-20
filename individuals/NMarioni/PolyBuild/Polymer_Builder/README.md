# Polymer Builder
**poly_runner.sh** - .sh script to run poly_build.py

**poly_build.py** - python script to create custom length/monomer random polymers
 - **NOTE**: header.txt/includes.txt should be edited to suite your .top needs
   - header.txt should path to desired forcefield
   - include.txt should include desired #include paths (ex: water model, ion model, etc)
 - **NOTE**: Line 19-22 rotate_check should be set to your needs
 - **NOTE**: Function MakeGro() creates the polymer.gro file by placing respecitive monomer units with some spacing. The spacing may need to be adjusted using delta_ar (see the comment for more info)
 - **NOTE**: MakeXXX() Functions have commented out #writefile() calls at the end. This can be useful for validation/bug fixing on small polymers

***IMPORTANT***: So far, this code has been tested with ethylene and acrylate backbones connecting the monomer units
                 Extra care should be taken to ensure the code works as intended if you do not have such simple backbone structures