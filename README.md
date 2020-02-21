# pyCRTM - python interface to CRTM.

## Bryan M. Karpowicz, Ph.D. - USRA/GESTAR/NASA 610.1 Global Modeling and Assimilation Office

This is a basic python interface to CRTM v2.3.0. I think it will generally serve my purposes and probably the purposes of other researchers needing a quick python accessible RT model.

The user interface is designed to be very similar to the python RTTOV interface. So, the user sets profiles, passes them to an object for a desired sensor, runs either the forward model/K-matrix, and pulls the brightness temperature/Jacobian/transmission/emissivity out of the object.  

The user interface shouldn't change all that much, and hopefully won't get too broken if/when features are added.  

This `README` has 4 parts:

1. Installation -- installing this library
2. Test/Examples -- describing test in testCases subdirectory
3. Python path etc -- how to use this library in a project.
4. Using the interface -- HOWTO/run through on how to use this interface

- Bryan Karpowicz -- August 16, 2019
---------------------------------------------------------------------------------------- 

## 1. Installation:
 
Done via the `setup_pycrtm.py` script  

Usage:
```
setup_pycrtm.py [-h] --install INSTALL --rtpath RTPATH --jproc JPROC [--arch ARCH] [--inplace]

setup_pycrtm.py: error: the following arguments are required: --install, --rtpath, --jproc 
```

### Required:
* `--install` -  Path where you want to install CRTM and coefficient files
*  `--rtpath`  -  Directory where the CRTM tarball from EMC ftp site will be cached/downloaded if not present. 
* `--jproc`   -  The number of threads to pass compiler

### Optional:
* `--arch` select compiler/environment gfortran (gcc) and ifort (intel) have been tested.
* `--inplace` this will skip the building of CRTM, but instead just compile pycrtm interface and link to CRTM library specified in `RTPATH`.

In addition to installing CRTM this script will patch the source fix some in in/out structures that cause gfortran to fail (gfortran.patch), and will re-organize
some null pointers to make the k-matrix thread safe for OpenMP (kmatrix.patch).  

Example to install CRTM in this directory under a subdirectory `lib/`:
```
./setup_pycrtm.py  --install $PWD/lib/ --rtpath $PWD/lib/ --jproc 1
```
Once completed:

* `$PWD/pycrtm.cpython-37m-PLATFORM.so` <-- (will always reside here) f2py interface compiled by setup 
* `$PWD/crtm.cfg`                       <-- Path where the CRTM coefficients are stored 

Following the example the CRTM will be installed here:

* `$PWD/lib/crtm/config.log`            <-- usual config associated with CRTM
* `$PWD/lib/crtm/include`               <-- path to all compiled fortran modules
* `$PWD/lib/crtm/lib/libcrtm.a`         <-- usual crtm static library
* `$PWD/lib/crtm/crtm_coef`             <-- path to crtm_coefficients
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
Right now w/ CRTM v2.3.0 there are differences between cases with zero cloud fraction, and with clouds turned off, so these test will fail, and generate plots.

## 3. Python path etc - needs work, but works for me at the moment: 

Right now things aren't setup to install into a user's Python path. What I typically plan on doing is set this as a git submodule, and bring it into a project and import the module using something like:
```Python
from pycrtm.pyCRTM import profilesCreate, pyCRTM
```
Or, do something like is done in the testCases scripts and insert the path, then import.
```Python
pycrtmDir = [this directory]
sys.path.insert(0,pycrtmDir)
from pyCRTM import pyCRTM, profilesCreate
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

