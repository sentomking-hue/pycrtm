#!/usr/bin/env python3
import os, h5py, sys, argparse 
import numpy as np
from matplotlib import pyplot as plt
from pyCRTM import pyCRTM, profilesCreate
 
def main(sensor_id,plotme,atten):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    # create 4 profiles for each of the 4 cases
    profiles = profilesCreate( 4, 92 )
    storedTb = []
    storedEmis = []
    # Make Reff large so correct number of streams is used
    Reff = 1000

    # populate the cases, and previously calculated Tb from crtm test program. 
       
    for i,c in enumerate(cases):
        h5 = h5py.File(os.path.join(thisDir,'data',c) , 'r')
        profiles.Angles[i,0] = 1.0
        profiles.Angles[i,1] = 0 
        profiles.Angles[i,2] = 100.0  # 100 degrees zenith below horizon.
        profiles.Angles[i,3] = 0.0 # zero solar azimuth 
        profiles.Angles[i,4] = 1.0 
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
        profiles.clouds[i,:,0,0] = cld 
        # note: For Moradi DDA cloud coefficients (only ones that work with active sensor)
        #       cloud effective radius is not used to lookup optical properties.
        #       HOWEVER, effective radius is used to determine the number of streams used
        #       by the radiative transfer. To avoid strange results, set the effective radius
        #       to a large size parameter to utilize the maximum number of streams.
 
        zzz = np.zeros(np.asarray(h5['cloudConcentration']).shape)
        cld = zzz
        cld[idx] = Reff 
        profiles.clouds[i,:,0,1] = cld 
        profiles.aerosols[i,:,0,0] = np.asarray(h5['aerosolConcentration'])
        profiles.aerosols[i,:,0,1] = np.asarray(h5['aerosolEffectiveRadius'])
        profiles.aerosolType[i] = h5['aerosolType'][()]
        profiles.cloudType[i] = 4
 
        zzz = np.zeros(np.asarray(h5['cloudConcentration']).shape)
        cld = zzz
        cld[idx]=1.0

        profiles.cloudFraction[i,:] = cld 
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
    crtmOb.output_aerosol_K = False
    crtmOb.output_cloud_K = True
    crtmOb.output_attenuated = atten
    crtmOb.loadInst()
    
    crtmOb.runK()
    Refl = crtmOb.ReflectivityAttenuated
    zz3=crtmOb.CloudEffectiveRadiusK
    zz4=crtmOb.CloudConcentrationK
    zz5=crtmOb.CloudFractionK
    Height = crtmOb.Height
    tol = 1e-6
    #! compute the mie parameter, 2.pi.reff/lambda
    # Reff in microns
    MieParameter = 2.0 * np.pi * Reff * crtmOb.wavenumber/10000.0

    #! determine the number of streams based on mie parameter as done in crtm
    if ( MieParameter < 0.01 ):
        nstreams = 2
    elif ( MieParameter < 1.0 ):
        nstreams = 4
    else:
        nstreams = 6

    if(nstreams == 6):
        print("nstreams based of reff is correct :",nstreams)
    else:
        print("nstreams based of reff is wrong :",nstreams)


    if(atten):
        tstMax3 = abs(zz3.max() - 0.0) < tol 
        tstMax4 = abs(zz4.max() - 0.0) < tol
        tstMax5 = abs(zz5.max() - 41.57776071153529) < tol
        tstMin3 = abs(zz3.min() - 0.0) < tol
        tstMin4 = abs(zz4.min() - -16.532941553040253) < tol
        tstMin5 = abs(zz5.min() - 0.0) < tol
        if(tstMax3 and tstMax4 and tstMin3 and tstMin4 and tstMin5):
            print('Yay! Min/Max reflectivities passed')
        else:
            print('Boo! something failed.')
            print('val, tol',abs(zz3.max() - 0), tol )
            print('val, tol',abs(zz4.max() - 41.57776071153529) , tol)
            print('val, tol',abs(zz5.max() - 0) , tol)
            print('val, tol',abs(zz3.min() - 0), tol )
            print('val, tol',abs(zz4.min() - -16.532941553040253) , tol)
            print('val, tol',abs(zz5.min() - 0.0) , tol)
    else:
        tstMax3 = abs(zz3.max() - 0.0) < tol 
        tstMax4 = abs(zz4.max() - 0.8685889638065037) < tol
        tstMax5 = abs(zz5.max() - 0.0) < tol
        tstMin3 = abs(zz3.min() - 0.0) < tol
        tstMin4 = abs(zz4.min() - 0.0) < tol
        tstMin5 = abs(zz5.min() - - 4.3429448190325175) < tol
        if(tstMax3 and tstMax4 and tstMin3 and tstMin4 and tstMin5):
            print('Yay! Min/Max reflectivities passed')
        else:
            print('Boo! something failed.')
            print('val, tol',abs(zz3.max() - 0), tol )
            print('val, tol',abs(zz4.max() - 0.8685889638065037) , tol)
            print('val, tol',abs(zz5.max() - 0) , tol)
            print('val, tol',abs(zz3.min() - 0), tol )
            print('val, tol',abs(zz4.min() - 0) , tol)
            print('val, tol',abs(zz5.min() - - 4.3429448190325175) , tol)
 
    if(plotme):
        for i,c in enumerate(cases):
            f,(ax_cld_refl0,ax_cld_refl1,ax_cld_refl2,ax_cld_conc) = plt.subplots( ncols=4,nrows=1,figsize=(24,5) )
            freq = crtmOb.frequencyGHz
            chans = []
            for ii,ff in enumerate(freq):
                iii=ii+1
                chans.append(iii)
            chans = np.asarray(chans)     
            idx, = np.where(profiles.clouds[i,:,0,0]>0)
            ax_cld_refl0.scatter(zz3[0,i,idx,0], Height[i,idx])
            ax_cld_refl1.scatter(zz4[0,i,idx,0],Height[i,idx])
            ax_cld_refl2.scatter(zz5[0,i,idx],Height[i,idx])
            ax_cld_conc.scatter(profiles.clouds[i,:,0,0],Height[i,:])
            ax_cld_refl0.set_ylabel('Height [km]')
            ax_cld_refl0.set_xlabel('Reflectivity Jacobian [dBz/$\mu$m]')
            ax_cld_refl1.set_xlabel('Reflectivity Jacobian [dBz/kg m$^{-2}$]')
            ax_cld_refl2.set_xlabel('Reflectivity Jacobian [dBz/Cloud Fraction]')
            ax_cld_refl0.set_ylim([0, Height.max()])
            ax_cld_refl1.set_ylim([0, Height.max()])
            ax_cld_refl2.set_ylim([0, Height.max()])
            ax_cld_conc.set_ylim([0, Height.max()])
            ax_cld_conc.set_xlabel('Concentration [kg m$^{-2}]$')   
            plt.tight_layout()
            print("note: Jacobian is zero for effective radius, because it isn't used in DDA scattering table.")
            plt.savefig(sensor_id+'_'+c+'_reflectivity_jacobian.png')

   

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = "Jacobian output test for CrIS NSR.")
    parser.add_argument('--plot',help="Plot Jacobians flag",dest='plotme',action='store_true')
    parser.add_argument('--attenuated',help="Use attenuated Reflectivity for Jacobian",dest='atten',action='store_true')
    a = parser.parse_args()
    sensor_id = 'cpr_cloudsat'
    main(sensor_id,a.plotme,a.atten)
 
