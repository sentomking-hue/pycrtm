import os, sys, glob, struct, configparser,netCDF4
import numpy as np 
def findLib(thisDir):
    """
    Find crtm library and check for shared library *.so. If so is found, return
    library path, otherwise return nothing to use static library.
    """
    cfg = configparser.ConfigParser()
    if ( os.path.exists( os.path.join(thisDir,'pyCRTM','pycrtm_setup.txt') ) ):
        pycrtm_setup_dir = os.path.join(thisDir,'pyCRTM','pycrtm_setup.txt')
    else:
        f = open(os.path.join(thisDir,'pyCRTM_JCSDA-2.0.1.dist-info','RECORD'))
        lines = f.readlines()
        for l in lines:
            if('pycrtm_setup.txt' in l):
                pycrtm_setup_dir = l.split('.txt')[0]
                pycrtm_setup_dir = pycrtm_setup_dir+'.txt'
    cfg.read( os.path.join(thisDir,pycrtm_setup_dir) )
    setupdir = cfg['Setup']['crtm_install']
    if( os.path.exists( os.path.join( setupdir, 'lib') ) ):
        libdir = os.path.join( setupdir, 'lib')
    elif( os.path.exists( os.path.join( setupdir, 'lib64') ) ):
        libdir = os.path.join( setupdir, 'lib64')
    else:
        libdir = os.path.join(setupdir,'') 
    so = glob.glob(os.path.join(libdir,'*.so'))
    if(len(so)>0):
        return libdir
    else:
        return ""
def findLibDyld(thisDir):
    """
    Find crtm library and check for shared library *.so. If so is found, return
    library path, otherwise return nothing to use static library.
    """
    cfg = configparser.ConfigParser()
    if ( os.path.exists( os.path.join(thisDir,'pyCRTM','pycrtm_setup.txt') ) ):
        pycrtm_setup_dir = os.path.join(thisDir,'pyCRTM','pycrtm_setup.txt')
    else:
        f = open(os.path.join(thisDir,'pyCRTM_JCSDA-2.0.1.dist-info','RECORD'))
        lines = f.readlines()
        for l in lines:
            if('pycrtm_setup.txt' in l):
                pycrtm_setup_dir = l.split('.txt')[0]
                pycrtm_setup_dir = pycrtm_setup_dir+'.txt'
    cfg.read( os.path.join(thisDir,pycrtm_setup_dir) )
    setupdir = cfg['Setup']['crtm_install']
    if( os.path.exists( os.path.join( setupdir, 'lib') ) ):
        libdir = os.path.join( setupdir, 'lib')
    elif( os.path.exists( os.path.join( setupdir, 'lib64') ) ):
        libdir = os.path.join( setupdir, 'lib64')
    else:
        libdir = os.path.join(setupdir,'') 
    so = glob.glob(os.path.join(libdir,'libcrtm.*'))
    if(len(so)>0):
        return libdir
    else:
        return ""



def setLD_LIBRARY_PATH(libdir):
    """
    For a given path add or create LD_LIBRARY_PATH and do os.execv thing to do it in the current
    python environment.
    """    
    #Set the LD_LIBRARY_PATH to make it possible to used shared object. 
    old_ld = os.environ.get("LD_LIBRARY_PATH")
    if old_ld:
        if(libdir not in os.environ["LD_LIBRARY_PATH"]):
            os.environ["LD_LIBRARY_PATH"] = old_ld + ":" + libdir
            os.execv(sys.argv[0], sys.argv)
    elif(len(libdir)>0):
        os.environ["LD_LIBRARY_PATH"] = libdir
        #os.execv(sys.argv[0], sys.argv)

def setDYLD_LIBRARY_PATH(libdir):
    """
    Warn user if dyld, that you need to set environment variable.
    """   
    #Set the LD_LIBRARY_PATH to make it possible to used shared object. 
    so = glob.glob(os.path.join(libdir,'libcrtm.*')) 
    old_ld = os.environ.get("DYLD_LIBRARY_PATH")
    if(len(so)>0):
        if old_ld:
            if(libdir not in os.environ["DYLD_LIBRARY_PATH"]):
                #os.environ["DYLD_LIBRARY_PATH"] = old_ld + ":" + libdir
                print('set DYLD_LIBRARY_PATH!')
                print('export DYLD_LIBRARY_PATH={}'.format(libdir))
                os.execv(sys.argv[0], sys.argv)
        elif(len(libdir)>0):
            #os.environ["DYLD_LIBRARY_PATH"] = libdir
             print('set DYLD_LIBRARY_PATH!')
             print('export DYLD_LIBRARY_PATH={}'.format(libdir))
             sys.exit()
             #os.execv(sys.argv[0], sys.argv)


def crtmLevelsToLayers( pLevels ):
    num = pLevels[1::] - pLevels[0:pLevels.shape[0]-1]
    den = np.log(pLevels[1::]/pLevels[0:pLevels.shape[0]-1])
    return num/den

def readTauCoeffODPS(fname):
    """
    Read and ODPS coefficient file.
    This code looks weird, because for whatever reason there's an extra 8 bytes after each binary record put into these things.
    So, you kind of have to go through the CRTM fortran and see look for each line things get read in, then advance by 8 bytes.
    input : file path to desired ODPS TauCoeff file
    output : dictionary of ODPS information (ODPS), dictionary of Optran information (Optran)
    """
    f = open(fname,'rb')
    o = {}
    
    # crtm binary header stuff
    version, magicNumber = struct.unpack('ii',f.read(struct.calcsize('ii')))
    pad = f.read(8)

    # ODPS file specific stuff 
    o['release'],o['version'] = struct.unpack('ii',f.read(struct.calcsize('ii')))
    f.read(8)

    o['algorithm'], = struct.unpack('i',f.read(struct.calcsize('i')))
    f.read(8)

    # dimensions for ODPS structure.
    o['n_Layers'], o['n_Components'], o['n_Absorbers'], o['n_Channels'], o['n_Coeffs'], o['n_OPIndex'], o['n_OCoeffs'] = struct.unpack('7i',f.read(struct.calcsize('7i')))
    f.read(8)
   
    n_Layers, n_Components, n_Absorbers, n_Channels, n_Coeffs, n_OPIndex, n_OCoeffs = o['n_Layers'], o['n_Components'], o['n_Absorbers'], o['n_Channels'], o['n_Coeffs'], o['n_OPIndex'], o['n_OCoeffs'] 
    # group index (not that useful)
    o['group_index'], = struct.unpack('i',f.read(struct.calcsize('i')))
    f.read(8)

    # sensor information.
    o['sensor_string'],o['wmo_satellite_id'], o['wmo_sensor_id'], o['sensor_type'] = struct.unpack('20s3i',f.read(struct.calcsize('20s3i')))
    f.read(8)

    # read in sensor channels in coef file
    fmt = '{:d}i'.format(n_Channels)
    o['sensor_channel'] = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    f.read(8)
    
    # read in components
    fmt = '{:d}i'.format(n_Components)
    o['component_id'] = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    f.read(8)

    # read in absorbers
    fmt = '{:d}i'.format(n_Absorbers)
    o['absorber_id'] = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    f.read(8)

    fmt = '{:d}d'.format(n_Layers+1)
    o['level_pressure'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    
    fmt = '{:d}d'.format(n_Layers)
    o['layer_pressure'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )

    fmt = '{:d}d'.format(n_Layers)
    o['layer_temperature'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )

    fmt = '{:d}d'.format(n_Layers*n_Absorbers)
    o['ref_absorber'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )


    fmt = '{:d}d'.format(n_Layers*n_Absorbers)
    o['min_absorber'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )

    fmt = '{:d}d'.format(n_Layers*n_Absorbers)
    o['max_absorber'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    f.read(8)

    # predictor and indexing stuff
    fmt = '{:d}i'.format(n_Channels*n_Components)
    o['n_predictors'] = struct.unpack(fmt, f.read( struct.calcsize(fmt) ) )

    fmt = '{:d}i'.format(n_Channels*n_Components)
    o['pos_index'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    f.read(8)
    
    # the actual ODPS coefficients     
    fmt = '{:d}f'.format(n_Coeffs)
    o['C'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    f.read(8)
    ODPS = o
    oo = {}
    # Old Optran coeff stuff.
    fmt = '{:d}i'.format(n_Channels)
    oo['OSignificance'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    order = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 

    fmt = '{:d}i'.format(n_Channels*7)
    oo['op_index'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    
    fmt = '{:d}i'.format(n_Channels)
    oo['op_pos_idx'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    
    fmt = '{:d}d'.format(n_OCoeffs)
    oo['OC'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    oo['alpha'], oo['alpha_c1'], oo['alpha_c2'], oo['oComponent_Index'] = struct.unpack('dddi',f.read(struct.calcsize('dddi')))
    oo['n_OCoeffs'] = n_OCoeffs
    Optran = oo
    f.close()

    return ODPS, Optran

def readNLTE(f,o):
    """
    Read non-local thermodynamic equilibrium coefficients (tied in with SpcCoeffs)
    input: file handle for spectral coefficient
    input/output: o containting input information, and output information
    """    

    o['release'],o['version'] = struct.unpack('ii',f.read(struct.calcsize('ii')))
    f.read(8)
    o['n_Predictors'], o['n_Sensor_Angles'] ,o['n_Solar_Angles']  , o['n_NLTE_Channels'] , o['n_Channels'] = struct.unpack('5i', f.read( struct.calcsize('5i') ))
    f.read(8)
    o['Sensor_Id'] = struct.unpack('20s', f.read( struct.calcsize('20s') ) )
    o['WMO_Satellite_Id'], o['WMO_Sensor_Id']  = struct.unpack('ii', f.read( struct.calcsize('ii'))) 
    fmt = '{:d}i'.format(o['n_Channels'])
    o['Sensor_Channel'] = struct.unpack(fmt, f.read( struct.calcsize(fmt) ) )
    f.read(8)
     
    fmt = '{:d}d'.format(2)
    o['Upper_Plevel'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    o['Lower_Plevel'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    f.read(8)

    o['Min_Tm'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    o['Max_Tm'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) )
    o['Mean_Tm'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    f.read(8)
    fmt = '{:d}i'.format(o['n_NLTE_Channels'])
    o['NLTE_Channel'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    f.read(8)
   
    fmt = '{:d}d'.format(o['n_Sensor_Angles'])
    o['Secant_Sensor_Zenith'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    fmt = '{:d}d'.format(o['n_Solar_Angles'])
    o['Secant_Solar_Zenith'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    f.read(8)
    
    fmt = '{:d}i'.format(o['n_Channels'])
    o['C_Index'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    fmt = '{:d}d'.format(o['n_Predictors']*o['n_Sensor_Angles']*o['n_Solar_Angles']*o['n_NLTE_Channels'])
    o['C'] = struct.unpack( fmt, f.read( struct.calcsize(fmt) ) ) 
    f.read(8)
    return o

def readSpcCoeff(fname):
    """
    Read Spectral Coefficient information.
    """
    f = open(fname,'rb')
    # crtm binary header stuff
    version, magicNumber = struct.unpack('ii',f.read(struct.calcsize('ii')))
    pad = f.read(8)
    o = {}
    # SpcCoeff file specific stuff 
    o['release'],o['version'] = struct.unpack('ii',f.read(struct.calcsize('ii')))
    f.read(8)
    o['n_Channels'], = struct.unpack('i',f.read(struct.calcsize('i')))
    n_Channels = o['n_Channels']
    f.read(8)
    # sensor information.
    o['sensor_string'], o['Sensor_Type'], o['wmo_satellite_id'], o['wmo_sensor_id'] = struct.unpack('20s3i',f.read(struct.calcsize('20s3i')))
    f.read(8)
    #information we probably care about.
    fmt = '{:d}i'.format(n_Channels)
    o['Sensor_Channel'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Polarization'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    if(o['version']>2):
        o['PolAngle'] = struct.unpack('{:d}d'.format(n_Channels),f.read(struct.calcsize('{:d}d'.format(n_Channels))))
    o['Channel_Flag'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))

    fmt = '{:d}d'.format(n_Channels)
    o['Frequency'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Wavenumber'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Planck_C1'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Planck_C2'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Band_C1'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Band_C2'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Cosmic_Background_Radiance'] = struct.unpack(fmt,f.read(struct.calcsize(fmt)))
    o['Solar_Irradiance'] =  struct.unpack(fmt,f.read(struct.calcsize(fmt))) 
    f.read(8)

    o['antenna_correction_present'], =  struct.unpack('i',f.read(struct.calcsize('i')))
    f.read(8)
    if(o['antenna_correction_present']==0): 
        o['nlte_correction_present'], = struct.unpack('i',f.read(struct.calcsize('i')))
        f.read(8)
        if(o['nlte_correction_present']>0): o = readNLTE(f,o)
    spcCoeff = o 
    f.close()

    return spcCoeff
def readSpcCoeffNc(fname):
    """
    Read Spectral Coefficient information from netcdf file.
    """
    o = {}
    ds = netCDF4.Dataset(fname,'r')
    for v in ds.variables:
        o[v] = np.asarray(ds.variables[v])

    for v in list(ds.ncattrs()):
        o[v.lower()] = np.asarray(ds.getncattr(v))
    return o

if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    # Stuff to get the installed rttov path, and import pyrttov interface
    pathInfo.read('crtm.cfg')
    spcCoeff = readSpcCoeff(os.path.join(pathInfo['CRTM']['coeffs_dir'],'cris399_npp.SpcCoeff.bin'))
    print('Spc Coeffs')
    for k in list(spcCoeff.keys()):
        print(k,spcCoeff[k])

    a, b = readTauCoeffODPS(os.path.join(pathInfo['CRTM']['coeffs_dir'],'cris399_npp.TauCoeff.bin'))
    print('ODPS')
    for k in list(a.keys()):
        print(k, a[k])
    print('OPTRAN')
    for k in list(b.keys()):
        print(k,b[k])
    
