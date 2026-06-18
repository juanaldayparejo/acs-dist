#! /bin/bash

###
### fmcd_gfortran.sh
###
### ---> This should be used to make MCD fortran stuff directly accessible in python 
### ---> This script is for gfortran, but it is easy to adapt to your other compilers
### ---> A file fmcd.so should be created
### ---> See mcd.py for use in python. Very easy!
###
### AS. 17/04/2012. 
### TP. 06/2022 adapted to MCD6.1

### PATHS to NETCDF and MCD software: Adapt these to your local settings
NETCDF="/usr"
wheremcd="/srv/workspace/data/mcd/MCD_6.1/mcd"
version="6.1"


### LOG FILE
\rm -f fmcd.log

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail (i.e. remove all '!' commments in codes)
sed s/"\!\!'"/"'"/g $wheremcd/MCD.F90            | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.MCD.F90

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
\rm -f fmcd.pyf
f2py -h fmcd.pyf -m fmcd tmp.MCD.F90 > fmcd.log 2>&1

#### BUILD
f2py -c fmcd.pyf -m fmcd tmp.MCD.F90 --fcompiler=gnu95 \
  -L$NETCDF/lib -lnetcdff -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-fPIC -ffree-form -ffree-line-length-none" \
  --verbose \
  > fmcd.log 2>&1

### CLEAN THE PLACE
\rm -f tmp.MCD.F90
