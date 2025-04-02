#!/usr/bin/env bash
module purge
# load top level module file to get conda (hopefully this doesn't disappear otherwise, we'll have to install that).
#module load miniforge
# use NCCS anaconda (not sure this is 100% neccessary, but works with jupyterhub)
module load anaconda
module use /gpfsm/dswdev/jcsda/spack-stack/scu17/spack-stack-1.9.0/envs/ue-intel-2021.10.0/install/modulefiles/Core
module load stack-intel/2021.10.0
module load stack-intel-oneapi-mpi/2021.10.0
module load cmake
module load git-lfs/3.4.0
export CONDA_VENV='pycrtm_intel'
conda create --name ${CONDA_VENV} python=3.11 ipykernel cartopy scipy dask xarray scikit-build h5py netcdf4 
conda activate ${CONDA_VENV}

export PWDOLD=${PWD}
export CHECKOUT_PATH="${PWD}/.local"
export PYCRTM_PATH="${PWD}"



#load modules here so have netcdf env variables in conda
module load netcdf-c/4.9.2 
module load netcdf-cxx4/4.3.1
module load netcdf-fortran/4.6.1 
# Activate the environment
# adjust this to your own tastes

# checkout and buildtest CRTM
mkdir $CHECKOUT_PATH
echo ${CHECKOUT_PATH}
cd ${CHECKOUT_PATH}
echo $PWD
git clone git@github.com:JCSDA/CRTMv3.git
cd CRTMv3
git checkout v3.1.1+build1
conda list
mkdir build 
cd build
cmake ../
make -j8 
make install 
ctest -j8 
#cd $CHECKOUT_PATH

cd ${PYCRTM_PATH}
#git clone git@github.com:JCSDA/pycrtm.git
#git checkout feature/prepCRTMv3
rm setup.cfg
printf "[Setup]\n">>setup.cfg
printf "# Specify the location of the crtm install\n">>setup.cfg 
printf "crtm_install =${CHECKOUT_PATH}/CRTMv3/build/\n" >> setup.cfg
printf "link_from_source_to_path_used = True\n" >> setup.cfg
printf "[Coefficients]\n" >> setup.cfg
printf "# source specify coefficient directory will grab little endian binary coefficients and netcdf and link them to path_used\n">>setup.cfg
printf "source_path =${CHECKOUT_PATH}/CRTMv3/build/test_data/\n" >> setup.cfg 
printf "# path used by pycrtm to read coefficients\n" >> setup.cfg 
printf "path_used =${CHECKOUT_PATH}/crtm_coefficients\n" >> setup.cfg 
python3 setup.py install
conda env config vars set LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
cd testCases
python3 test_atms.py 
cd $PWDOLD 
