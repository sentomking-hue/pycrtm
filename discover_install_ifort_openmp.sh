#!/bin/sh
source $PWD/discover_modules_ifort_openmp.sh 
./setup_pycrtm.py --coef $PWD --install $PWD/../REL-2.4.0-alpha/  --repos $PWD/../REL-2.4.0-alpha/ --jproc 2  --ncpath $NETCDF_ROOT --h5path $HDF5_ROOT --arch ifort-openmp
