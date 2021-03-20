# pyCRTM - python interface to CRTM.

## Bryan M. Karpowicz, Ph.D. - USRA/GESTAR/NASA 610.1 Global Modeling and Assimilation Office, with Contributions from Patrick Stegmann, Dr.-Ing. - JCSDA

This is a basic python interface to CRTM v2.4.0. 

The user interface is designed to be very similar to the python RTTOV interface. So, the user sets profiles, passes them to an object for a desired sensor, runs either the forward model/K-matrix, and pulls the brightness temperature/Jacobian/transmission/emissivity out of the object.  


This `README` has 4 parts:

1. Installation -- installing this library
2. Test/Examples -- describing test in testCases subdirectory
3. Importing -- how to use this library in a project.
4. Using the interface -- HOWTO/run through on how to use this interface

- Bryan Karpowicz -- October 23, 2020
---------------------------------------------------------------------------------------- 

## 1. Installation:
A. Dependencies numpy and scikit-build (install those first, if you don't have them.) 
B. Done via the `setup.py` invoking standard setuptools-style 
```
python3 setup.py install
```
compiler options are handled through setup.cfg
```
[Setup]
# Options
# compiler = gfortran
# compiler = intel
# compiler = gfortran-openmp
# compiler = intel-openmp
compiler = gfortran-openmp
crtm_install = /discover/nobackup/bkarpowi/github/JCSDA_crtm/crtm-bundle/crtm_v2.4.0/
# Download Coefficients
# Controls whether coefficients are downloaded
download = True
coef_with_install = True
[Coefficients]
# Use to specify alternative coefficient file location where Little Endian Coefficient files are stored.
# If user desires coefficients to be stored with the installed package, leave this alone.
# If user selects coef_with_install = False, this must be specified.
# set argument below (path) to the full path of the coefficients.
path = /discover/nobackup/projects/gmao/obsdev/bkarpowi/tstCoef/
```
By default the CRTM coefficients will be copied along with the pycrtm install to your site-packages. To change this set `coef_with_install` to False, and set `path` to the desired location. Note: `path` will beignored if `coef_with_install` is True.

If you don't have write access to your python distribution you can use the standard --prefix or --user to install to a local directory paired with setting your `PYTHONPATH` to point to the local directory.
 

---------------------------------------------------------------------------------------- 

## 2. Tests/Examples:

A few test cases have been developed using input/output grabbed from the CRTM test program.
Two basic scripts which will perform cases 1-4 (stored in $PWD/testCases/data/case[n].h5) from the CRTM test program on 4 OpenMP threads: 
* `$PWD/testCases/test_atms.py`
* `$PWD/testCases/test_cris.py`
These *should* just say Yay, and not produce any plots if successful. 

The following scripts will do the same thing, only this time load up the same 4 profiles multiple times to further test threading with 10 threads (turn your laptop into a space heater more or less).
* `$PWD/testCases/test_atms_threads.py`
* `$PWD/testCases/test_cris_threads.py`
These *should* just say Yay, and not produce any plots if successful. 


The following scripts will run CRTM without aerosols or clouds:
* `$PWD/testCases/test_atms_no_clouds.py`
* `$PWD/testCases/test_cris_no_clouds.py`

For those Jupyter notebook fans, there is even Jupyter notebook example simulating ATMS:
* `$PWD/testCases/test_atms.ipynb`

## 3. Importing 

```Python
from pycrtm.pyCRTM import profilesCreate, pyCRTM
```
---------------------------------------------------------------------------------------- 

## 4. Using the interface (designed to be pretty much like the RTTOV equivalent python library):

Create a profiles data structure using `profilesCreate(nprofiles, nlayers)` which will generate an object with user specified number of profiles, and number of layers in the profiles provided.
```Python
profiles = profilesCreate(4, 92) # will generate an empty object with 4 profiles each with 92 layers. 
```
Once initialized, the user will need to provide values for the desired profiles (see example scripts). Next, the user initializes a crtm instance, set desired parameters, and passes profiles to the CRTM:

```Python
crtmOb = pyCRTM()
crtmOb.coefficientPath = pathInfo['CRTM']['coeffs_dir']
crtmOb.sensor_id = sensor_id
crtmOb.nThreads = 4
crtmOb.profiles = profiles
```

Next, the instrument is loaded/number of channels in the output structure are initialized:

```Python
crtmOb.loadInst()
```

Next, the user can run either the forward model (runDirect), or the K-matrix (runK) 
```Python
crtmOb.runDirect()
crtmOb.runK()
```
Finally, the user can pull out desired parameters such as brightness temperatures, Jacobians, or Transmission along path (TauLevels - the derivative will be the weighting function):
```Python
# brightness temperature (nprofiles, nchan):
brightnessTemperature = crtmOb.Bt 

#Transmission (to compute weighting functions) ( nprofiles, nchan, nlayers)
Tau = crtmOb.TauLevels 

#Temperature, Water Vapo[u]r, and Ozone Jacobians ( npforfiles, nchan, nlayers)
O3_Jacobian = crtmOb.O3K
Water_Vapor_Jacobian = crtmOb.QK
Temperature_Jacobian = crtm.TK

#Emissivity (nprofiles, nchan)
Emissivity = crtmOb.surfEmisRefl
```
---------------------------------------------------------------------------------------- 

