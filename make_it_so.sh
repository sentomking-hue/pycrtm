#!/bin/sh
export SEL=$1

if [[ ${SEL} == 'apple_silicon' ]]; then
    mkdir -p ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
    zsh ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 
    rm ~/miniconda3/miniconda.sh
    source ~/miniconda3/bin/activate
    conda init --all
elif [[ ${SEL} == 'apple_intel' ]]; then
    mkdir -p ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    source ~/miniconda3/bin/activate
    conda init --all
elif [[ ${SEL} == 'linux' ]]; then
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    source ~/miniconda3/bin/activate
    conda init --all
elif [[ ${SEL} == 'skip' ]]; then
    echo "Skipping miniconda install. Already installed."    
else
    echo "Unknown platform selected:${1}"
	exit 1
fi
export CONDA_VENV='pycrtm'
conda create --name ${CONDA_VENV} python=3.11 scikit-build h5py netcdf4 gfortran libnetcdf netcdf-fortran cmake git git-lfs matplotlib
conda init
conda activate ${CONDA_VENV}

export PWDOLD=${PWD}
export CHECKOUT_PATH="${PWD}/ext"
export PYCRTM_PATH="${PWD}"


# checkout and buildtest CRTM
mkdir $CHECKOUT_PATH
echo ${CHECKOUT_PATH}
cd ${CHECKOUT_PATH}
echo $PWD
git clone https://github.com/JCSDA/CRTMv3.git
cd CRTMv3
git checkout v3.1.1+build1
conda list
mkdir build 
cd build
cmake ../
make -j8 
make install 
ctest -j8 

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

cd testCases
if [[ $SEL == 'apple_intel' ]]; then
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${CHECKOUT_PATH}/CRTMv3/build/lib/" 
elif [[ $SEL == 'skip' ]]; then
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${CHECKOUT_PATH}/CRTMv3/build/lib/" 
fi
python3 test_atms.py 
cd $PWDOLD 
