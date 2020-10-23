#!/bin/sh
module purge
set OPT=/discover/swdev/jcsda/modules
module use $OPT/modulefiles/apps
module use $OPT/modulefiles/core
module load jedi/gnu-impi/9.2.0
module load python/GEOSpyD/Ana2019.10_py3.7
