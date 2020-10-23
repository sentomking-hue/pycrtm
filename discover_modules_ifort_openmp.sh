#!/bin/sh
module purge
set OPT=/discover/swdev/jcsda/modules
module use $OPT/modulefiles/apps
module use $OPT/modulefiles/core
module load jedi/intel-impi/19.1.0.166
module load python/GEOSpyD/Ana2019.10_py3.7
