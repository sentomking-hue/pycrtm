module purge
module use /discover/swdev/jcsda/spack-stack/modulefiles
module load miniconda
module load ecflow
module load mysql
module use /gpfsm/dswdev/jcsda/spack-stack/spack-stack-1.5.0/envs/unified-env/install/modulefiles/Core
module load stack-gcc
module load stack-openmpi
module load stack-python
module load jedi-fv3-env
module load py-matplotlib
module load py-cartopy
module list
# add your python path here: e.g.
#setenv PYTHONPATH ${PYTHONPATH}:/discover/nobackup/projects/gmao/obsdev/bkarpowi/pycrtm_builds/pythonModules3cmake
