#!/usr/bin/env python3
import os, h5py, sys, argparse 
import numpy as np
from matplotlib import pyplot as plt
from pyCRTM import pyCRTM, profilesCreate
 
def main(sensor_id,plotMe):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    # create 4 profiles for each of the 4 cases
    profiles = profilesCreate( 4, 92)
    storedTb = []
    storedEmis = []
    # populate the cases, and previously calculated Tb from crtm test program.    
    for i,c in enumerate(cases):
        h5 = h5py.File(os.path.join(thisDir,'data',c) , 'r')
        profiles.Angles[i,0] = 1.0#h5['zenithAngle'][()]
        profiles.Angles[i,1] = 0 
        profiles.Angles[i,2] = 100.0  # 100 degrees zenith below horizon.
        profiles.Angles[i,3] = 0.0 # zero solar azimuth 
        profiles.Angles[i,4] = 1.0 #h5['scanAngle'][()]
        profiles.DateTimes[i,0] = 2001
        profiles.DateTimes[i,1] = 1
        profiles.DateTimes[i,2] = 1
        profiles.Pi[i,:] = np.asarray(h5['pressureLevels'] )
        profiles.P[i,:] = np.asarray(h5['pressureLayers'][()])
        profiles.T[i,:] = np.asarray(h5['temperatureLayers'])
        profiles.Q[i,:] = np.asarray(h5['humidityLayers'])
        profiles.O3[i,:] = np.asarray(h5['ozoneConcLayers'])
        idx = np.where(np.asarray(h5['cloudConcentration'])>0)
        zzz = np.zeros(np.asarray(h5['cloudConcentration']).shape)
        cld = zzz
        cld[idx] = 5
        profiles.clouds[i,:,0,0] = cld #np.asarray(h5['cloudConcentration'])
        zzz = np.zeros(np.asarray(h5['cloudConcentration']).shape)
        cld = zzz
        cld[idx] = 1000 
        profiles.clouds[i,:,0,1] = cld #np.asarray(h5['cloudEffectiveRadius'])
        profiles.aerosols[i,:,0,0] = np.asarray(h5['aerosolConcentration'])
        profiles.aerosols[i,:,0,1] = np.asarray(h5['aerosolEffectiveRadius'])
        profiles.aerosolType[i] = h5['aerosolType'][()]
        profiles.cloudType[i] = 4# h5['cloudType'][()]
 
        zzz = np.zeros(np.asarray(h5['cloudConcentration']).shape)
        cld = zzz
        cld[idx]=1.0

        profiles.cloudFraction[i,:] = cld #h5['cloudFraction'][()]
        profiles.climatology[i] = h5['climatology'][()]
        profiles.surfaceFractions[i,:] = h5['surfaceFractions']
        profiles.surfaceTemperatures[i,:] = h5['surfaceTemperatures']
        profiles.Salinity[i] = 33.0 
        profiles.windSpeed10m[i] = 5.0
        profiles.LAI[i] = h5['LAI'][()]
        profiles.windDirection10m[i] = h5['windDirection10m'][()]
        # land, soil, veg, water, snow, ice
        profiles.surfaceTypes[i,0] = h5['landType'][()]
        profiles.surfaceTypes[i,1] = h5['soilType'][()]
        profiles.surfaceTypes[i,2] = h5['vegType'][()]
        profiles.surfaceTypes[i,3] = h5['waterType'][()]
        profiles.surfaceTypes[i,4] = h5['snowType'][()]
        profiles.surfaceTypes[i,5] = h5['iceType'][()]
        storedTb.append(np.asarray(h5['Tb_atms'][0:22]))
        storedEmis.append(np.asarray(h5['emissivity_atms'][0:22]))
        h5.close()

    crtmOb = pyCRTM()
    crtmOb.profiles = profiles
    crtmOb.sensor_id = sensor_id
    crtmOb.nThreads = 1
    crtmOb.Active = True
    crtmOb.CloudCoeff_File = 'CloudCoeff_DDA_Moradi_2022.nc4'
    crtmOb.loadInst()
    
    crtmOb.runDirect()
    forwardReflectivity = crtmOb.Reflectivity
    forwardReflectivityAttenuated = crtmOb.ReflectivityAttenuated
    Height= crtmOb.Height
    zz3 = forwardReflectivity
    zz4 = forwardReflectivityAttenuated 
    idx3 = np.where(zz3>-9000)
    idx4 = np.where(zz4>-9000)
    tol = 1e-6
    tstMax3 = abs(zz3[idx3].max() - 30.21725448020621) < tol 
    tstMax4 = abs(zz4[idx4].max() - 24.95743632633681) < tol
    tstMin3 = abs(zz3[idx3].min() - 29.50054261447972) < tol
    tstMin4 = abs(zz4[idx4].min() - - 59.492652913112856) < tol
    if(tstMax3 and tstMax4 and tstMin3 and tstMin4):
        print('Yay! Min/Max reflectivities passed')
    else:
        print('Boo! something failed.')
        print('val, tol',abs(zz3[idx3].max() - 30.21725448020621), tol )
        print('val, tol',abs(zz4[idx4].max() - 24.95743632633681) , tol)
        print('val, tol',abs(zz3[idx3].min() - 29.50054261447972) , tol)
        print('val, tol',abs(zz4[idx4].min() - - 59.492652913112856) , tol)
    if(plotMe):
        for i,c in enumerate(cases):
            f,(ax_cld_refl0,ax_cld_refl1,ax_cld_conc) = plt.subplots( ncols=3,nrows=1,figsize=(24,5) )
            freq = crtmOb.frequencyGHz
            chans = []
            for ii,ff in enumerate(freq):
                iii=ii+1
                chans.append(iii)
            chans = np.asarray(chans)     
            pgrid,nu_grid = np.meshgrid(chans,profiles.P[i,:])
            idx, = np.where(zz3[i,0,:]>-9000)
            ax_cld_refl0.scatter(zz3[i,0,idx], Height[i,idx])
            idx, = np.where(zz4[i,0,:]>-9000)
            ax_cld_refl1.scatter(zz4[i,0,idx],Height[i,idx])
            ax_cld_conc.scatter(profiles.clouds[i,:,0,0],Height[i,:])
            ax_cld_refl0.set_ylabel('Height [km]')
            ax_cld_refl0.set_xlabel('Reflectivity [dBz]')
            ax_cld_refl1.set_xlabel('Reflectivity Attenuated [dBz]')
            ax_cld_refl0.set_ylim([0, Height.max()])
            ax_cld_refl1.set_ylim([0, Height.max()])
            ax_cld_conc.set_ylim([0, Height.max()])
            ax_cld_conc.set_xlabel('Concentration [kg m$^{-2}]$')   
            plt.tight_layout()
            plt.savefig(sensor_id+'_'+c+'_reflectivity.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = "Jacobian output test for CrIS NSR.")
    parser.add_argument('--plot',help="Plot Jacobians flag",dest='plotme',action='store_true')
    a = parser.parse_args()
    sensor_id = 'cpr_cloudsat'
    main(sensor_id,a.plotme)
 
