#!/usr/bin/env python3
import os, h5py, sys, argparse 
import numpy as np
from matplotlib import pyplot as plt
from pyCRTM import pyCRTM, profilesCreate
 
def main(sensor_id,plotme):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    # create 4 profiles for each of the 4 cases
    profiles = profilesCreate( 4, 92 )
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
        cld[idx]=1000 
        profiles.clouds[i,:,0,1] =cld #np.asarray(h5['cloudEffectiveRadius'])
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
    crtmOb.output_aerosol_K = False
    crtmOb.output_cloud_K = True
    crtmOb.loadInst()

    
    crtmOb.runK()
    Refl = crtmOb.ReflectivityAttenuated
    print('whir0')
    print(Refl)
    zz3=crtmOb.CloudEffectiveRadiusK
    zz4=crtmOb.CloudConcentrationK
    zz5=crtmOb.CloudFractionK
    print('whir1')
    print(zz3)
    print("whir2")
    print(zz4)

    print("whir3")
    print(zz5)
    """
    if(tst3 and tst4 and tst5):
        print("Yay! Cloud Jacobians Pass.")
    else:
       print("Boo! Cloud Jacobians Fail.")
       print("Effective Radius Test",tst3)
       print("Concentratoin Test",tst4)
       print("Fraction Test",tst5)
    """
    if(plotme):
        for i,c in enumerate(cases):
            #ignore first two profiles, because cloud fraction set to zero.
            if(i>1):
                f,(ax_cld_jac0,ax_cld_jac1,ax_cld_jac2,ax_cld_conc) = plt.subplots( ncols=4,nrows=1,figsize=(24,5) )
                chans = np.arange(1,23,1)
                freq = crtmOb.frequencyGHz
                pgrid,nu_grid = np.meshgrid(chans,profiles.P[i,:])
                s1 = ax_cld_jac0.scatter(pgrid,nu_grid,c=zz3[:,i,:,0].T,s=3)
                s2 = ax_cld_jac1.scatter(pgrid,nu_grid,c=zz4[:,i,:,0].T,s=3 )
                s3 = ax_cld_jac2.scatter(pgrid,nu_grid,c=zz5[:,i,:].T,s=3 )
                ax_cld_conc.scatter(profiles.clouds[i,:,0,0],profiles.P[i,:])
                #ax_cld_conc.scatter(profiles.cloudFraction[i,:],profiles.P[i,:])

                ax_cld_conc.invert_yaxis()
                ax_cld_jac0.invert_yaxis()
                ax_cld_jac1.invert_yaxis()
                ax_cld_jac2.invert_yaxis()
                ax_cld_jac0.set_ylabel('Pressure [hPa]')
                ax_cld_jac0.set_xlabel('Channel')
                ax_cld_jac1.set_xlabel('Channel')
                ax_cld_jac2.set_xlabel('Channel')

                ax_cld_jac0.set_xticks(chans)
                ax_cld_jac1.set_xticks(chans)
                ax_cld_jac2.set_xticks(chans)
                ax_cld_conc.set_xlabel('Concentration [kg m$^{-2}]$')   
                f.colorbar(s1,label='R$_e$ [$\mu$m K$^{-1}$]')
                f.colorbar(s2,label='Concentration Jacobian [kg m$^{-2}$ K$^{-1}$]')
                f.colorbar(s3,label='Fraction Jacobian [K$^{-1}$]')
                plt.tight_layout()
                plt.savefig(sensor_id+'_'+c+'_cloud_jacobian.png')
       
  

    #if ( all( np.abs( forwardTb.flatten() - np.asarray(storedTb).flatten() ) <= 1e-5)  and all( np.abs( kTb.flatten() - np.asarray(storedTb).flatten() ) <= 1e-5) ):
    #    print("Yay! all values are close enough to what CRTM test program produced!")
    #else: 
    #    print("Boo! something failed. Look at cris plots")
    #    wavenumbers = np.zeros([4,1305])
    #    wavenumbers[0:4,:] = np.linspace(1,1306,1305)
    #    plt.figure()
    #    plt.plot(wavenumbers.T,forwardTb.T-np.asarray(storedTb).T ) 
    #    plt.legend(['1','2','3','4'])
    #    plt.savefig(os.path.join(thisDir,'cris'+'_spectrum_forward.png'))
    #    plt.figure()
    #    plt.plot(wavenumbers.T,forwardEmissivity.T-np.asarray(storedEmis).T)
    #    plt.savefig(os.path.join(thisDir,'cris'+'_emissivity_forward.png')) 
    #
    #    plt.figure()
    #    plt.plot(wavenumbers.T,kTb.T-np.asarray(storedTb).T)
    #    plt.savefig(os.path.join(thisDir,'cris'+'_spectrum_k.png'))
    #    plt.figure()
    #    plt.plot(wavenumbers.T,kEmissivity.T-np.asarray(storedEmis).T)
    #    plt.savefig(os.path.join(thisDir,'cris'+'_emissivity_k.png')) 
    #    sys.exit("Boo! didn't pass tolerance with CRTM test program.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = "Jacobian output test for CrIS NSR.")
    parser.add_argument('--plot',help="Plot Jacobians flag",dest='plotme',action='store_true')
    a = parser.parse_args()
    sensor_id = 'cpr_cloudsat'
    main(sensor_id,a.plotme)
 
