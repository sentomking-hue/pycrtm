# -*- makefile -*
#
#  Trying to run make manually? You'll need to set the following environment variables.
#  FORT -- Compiler you want (ifort, gfortran, etc)
#  F2PY_COMPILER -- what f2py calls the compiler (gfortran = gnu95, intel = intelem)
#  FCFLAGS -- flags taken for the appropriate compiler (see CRTM "config-setup" directory)
#  ILOC -- path to the CRTM install directory. 
#          The install directory where you find  "lib" (where the libcrtm.a lives) and "include" (where all the *.mod files live)
  F2PY = f2py --fcompiler=${F2PY_COMPILER} --f90flags='${FCFLAGS}'
  LIBA = -L${ILOC}/lib #-L /Users/bkarpowi/homebrew/Cellar/netcdf/4.7.4_1/lib -L /Users/bkarpowi/homebrew/Cellar/hdf5/1.12.0_1/lib 
  INCA = -I${ILOC}/include# -I /Users/bkarpowi/homebrew/Cellar/netcdf/4.7.4_1/include -I /Users/bkarpowi/homebrew/Cellar/hdf5/1.12.0_1/include
MODULE=pycrtm

all: ${MODULE}.so

#Only really need first bit, if you change interface, but do it anyway so you don't forget.
${MODULE}.so: pycrtm.f90
	f2py -m ${MODULE} -h sgnFile.pyf pycrtm.f90 --overwrite-signature
	${F2PY} ${F2PY_FLAGS}-c ${LIBA} -lcrtm -lgomp -L/Users/bkarpowi/homebrew/Cellar/netcdf/4.7.4_1/lib -lnetcdf -lnetcdff  -L/Users/bkarpowi/homebrew/Cellar/hdf5/1.12.0_1/lib  -lhdf5 ${INCA} -I/Users/bkarpowi/homebrew/Cellar/netcdf/4.7.4_1/include -I/Users/bkarpowi/homebrew/Cellar/hdf5/1.12.0_1/include -m ${MODULE} $<  only: wrap_forward wrap_k_matrix

clean:
	${RM} ${MODULE}*.so
