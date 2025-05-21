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
    #mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    source ~/miniconda3/bin/activate
    conda init --all
elif [[ ${SEL} == 'skip' ]]; then
    echo "Skipping miniconda install. Already installed."    
    if [[ -f "$HOME/miniconda3/bin/activate" ]]; then
        source $HOME/miniconda3/bin/activate
    fi
else
    echo "Unknown platform selected:${1}"
	exit 1
fi

# Add conda-forge, necessary to get netcdf-fortran
#conda config --add channels conda-forge

export CONDA_VENV='pycrtm'
conda create --name ${CONDA_VENV} python scikit-build h5py netcdf4 gfortran libnetcdf netcdf-fortran meson cmake git git-lfs matplotlib
conda activate ${CONDA_VENV}
conda init
echo $CONDA_PREFIX



export PWDOLD=${PWD}
export CHECKOUT_PATH="${PWD}/ext"
export PYCRTM_PATH="${PWD}"
export CRTM_VERSION="3.x"

# checkout and buildtest CRTM
mkdir ${CHECKOUT_PATH}
echo ${CHECKOUT_PATH}
cd ${CHECKOUT_PATH}
echo ${PWD}

if [[ ${CRTM_VERSION} == "2.x" ]]; then
  git clone https://github.com/ecmwf/ecbuild
  export PATH=${PWD}/ecbuild/bin:$PATH
  git clone https://github.com/JCSDA/crtm/
  cd crtm
  git checkout v2.4.1-jedi.1
fi

if [[ ${CRTM_VERSION} == "3.x" ]]; then
  git clone https://github.com/JCSDA/CRTMv3.git
  cd CRTMv3
  git checkout v3.1.1+build1
fi

export FPIC="-fPIC"
# make sure -fPIC is set before CRTM is built, add to FFLAGS
if [[ -z "$FFLAGS" ]]; then
  export FFLAGS="${FPIC}"
  echo "Environment variable FFLAGS not found. Setting it to '${FPIC}'."
elif [[ ! "$!FFLAGS" =~ (^|[[:space:]])"${FPIC}"($|[[:space:]]) ]]; then
  export FFLAGS="${FFLAGS} ${FPIC}"
  echo "Flag '${FPIC}' not found in FFLAGS. Appending it."
else
  echo "Flag '${FPIC}' is already set in FFLAGS."
fi

#make sure -fPIC is set before CRTM is build, add to FORTRANFLAGS
if [[ -z "$FORTRANFLAGS" ]]; then
  export FORTRANFLAGS="${FPIC}"
  echo "Environment variable FORTRANFLAGS not found. Setting it to '${FPIC}'."
elif [[ ! "$!FORTRANFLAGS" =~ (^|[[:space:]])"${FPIC}"($|[[:space:]]) ]]; then
  export FORTRANFLAGS="${FORTRANFLAGS} ${FPIC}"
  echo "Flag '${FPIC}' not found in FORTRANFLAGS. Appending it."
else
  echo "Flag '${FPIC}' is already set in FORTRANFLAGS."
fi
#conda list
if [[ ${CRTM_VERSION} == "2.x" ]]; then
  chmod +x Get_CRTM_Binary_Files.sh
  ./Get_CRTM_Binary_Files.sh
fi

mkdir build 
cd build
cmake -DBUILD_SHARED_LIBS=OFF ../
make -j8 
make install 
ctest -j8 

cd ${PYCRTM_PATH}
#git clone git@github.com:JCSDA/pycrtm.git
#git checkout feature/prepCRTMv3
rm setup.cfg
printf "[Setup]\n">>setup.cfg
printf "# Specify the location of the crtm install\n">>setup.cfg 

if [[ ${CRTM_VERSION} == "2.x" ]]; then
  printf "crtm_install =${CHECKOUT_PATH}/crtm/build/\n" >> setup.cfg
fi
if [[ ${CRTM_VERSION} == "3.x" ]]; then
  printf "crtm_install =${CHECKOUT_PATH}/CRTMv3/build/\n" >> setup.cfg
fi

printf "link_from_source_to_path_used = True\n" >> setup.cfg
printf "[Coefficients]\n" >> setup.cfg
printf "# source specify coefficient directory will grab little endian binary coefficients and netcdf and link them to path_used\n">>setup.cfg

if [[ ${CRTM_VERSION} == "2.x" ]]; then
  printf "source_path =${CHECKOUT_PATH}/crtm/fix/\n" >> setup.cfg 
fi

if [[ ${CRTM_VERSION} == "3.x" ]]; then
  printf "source_path =${CHECKOUT_PATH}/CRTMv3/build/test_data/\n" >> setup.cfg 
fi

printf "# path used by pycrtm to read coefficients\n" >> setup.cfg 
printf "path_used =${CHECKOUT_PATH}/crtm_coefficients\n" >> setup.cfg 

pip install . -v

cd testCases
if [[ $SEL == 'apple_intel' ]]; then
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${CHECKOUT_PATH}/CRTMv3/build/lib/" 
elif [[ $SEL == 'skip' ]]; then
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${CHECKOUT_PATH}/CRTMv3/build/lib/" 
fi
python3 test_atms.py 
cd $PWDOLD 
