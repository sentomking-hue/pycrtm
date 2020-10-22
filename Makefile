# -*- makefile -*
#
#  Trying to run make manually? You'll need to set the following environment variables.
#  FORT -- Compiler you want (ifort, gfortran, etc)
#  F2PY_COMPILER -- what f2py calls the compiler (gfortran = gnu95, intel = intelem)
#  FCFLAGS -- flags taken for the appropriate compiler (see CRTM "config-setup" directory)
#  ILOC -- path to the CRTM install directory. 
#          The install directory where you find  "lib" (where the libcrtm.a lives) and "include" (where all the *.mod files live)
#  NCINSTALL -- path to NETCDF install directory (above 'lib' and 'include')
#  H5INSTALL -- path to HDF5 install directory (abvoe 'lib' and 'include')
#  DASHL -- things like -lcrtm -lgomp -lnetcdf -lnetcdff, etc
  F2PY = f2py --fcompiler=${F2PY_COMPILER} --f90flags='${FCFLAGS}'
  LIBCRTM= -L${ILOC}/lib 
  INCCRTM = -I${ILOC}/include
  LIBNC = -L${NCINSTALL}/lib
  INCNC = -L${NCINSTALL}/include
  LIBH5 = -L${H5INSTALL}/lib
  INCH5 = -L${H5INSTALL}/include
MODULE=pycrtm

all: ${MODULE}.so

#Only really need first bit, if you change interface, but do it anyway so you don't forget.
${MODULE}.so: pycrtm.f90
	f2py -m ${MODULE} -h sgnFile.pyf pycrtm.f90 --overwrite-signature
	${F2PY} ${F2PY_FLAGS}-c ${LIBCRTM} ${LIBNC} ${LIBH5} ${INCCRTM} ${INCNC} ${INCH5} ${DASHL} -m ${MODULE} $<  only: wrap_forward wrap_k_matrix

clean:
	${RM} ${MODULE}*.so
