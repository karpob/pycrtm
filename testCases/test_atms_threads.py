#!/usr/bin/env python3
import configparser
import os, h5py, sys 
import numpy as np
from matplotlib import pyplot as plt
thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(thisDir)
sys.path.insert(0,parentDir)
from pycrtm import pycrtm
import time 
def main(coefficientPath, sensor_id):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    salinity = 33.0
    zenithAngle = []
    azimuthAngle = []
    scanAngle = []
    solarAngle = []
    pressureLevels = []
    pressureLayers = []
    temperatureLayers =[]
    humidityLayers = []
    ozoneConcLayers =[]
    co2ConcLayers = []
    aerosolEffectiveRadius = []
    aerosolConcentration = []
    aerosolType = []
    cloudEffectiveRadius = []
    cloudConcentration = []
    cloudType = []
    cloudFraction =[]
    climatology = []
    surfaceTemperatures = []
    surfaceFractions = []
    LAI = []
    windSpeed = []
    windDirection = []
    n_absorbers = []
    landType = []
    soilType = []
    vegType = []
    waterType = []
    snowType = []
    iceType = []
    duplicateCases = []
    storedTb = []
    storedEmis = []
    year = []
    month = []
    day = []
    ii = 0
    for i in list(range(200)):
        duplicateCases.append(cases[ii])
        ii+=1
        if(ii==4): ii=0
 
    for c in duplicateCases:
        #print(pycrtm.__doc__) 
        h5 = h5py.File(os.path.join(thisDir, 'data',c) , 'r')
        nChan = 22
        zenithAngle.append(h5['zenithAngle'][()])
        scanAngle.append(h5['scanAngle'][()])
        azimuthAngle.append( 999.9), 
        solarAngle.append(np.asarray([100.0,0.0]))
        year.append(2001)
        month.append(1)
        day.append(1)
        pressureLevels.append(h5['pressureLevels'])
        pressureLayers.append(h5['pressureLayers'])
        temperatureLayers.append(h5['temperatureLayers']) 
        humidityLayers.append(h5['humidityLayers']) 
        ozoneConcLayers.append(h5['ozoneConcLayers'])
        co2ConcLayers.append(h5['co2ConcLayers'])
        aerosolEffectiveRadius.append(h5['aerosolEffectiveRadius'])
        aerosolConcentration.append( h5['aerosolConcentration'] ) 
        aerosolType.append(h5['aerosolType'][()])
        cloudEffectiveRadius.append(h5['cloudEffectiveRadius']) 
        cloudConcentration.append(h5['cloudConcentration']) 
        cloudType.append(h5['cloudType'][()] ) 
        cloudFraction.append(h5['cloudFraction'])  
        climatology.append(h5['climatology'][()])
        surfaceTemperatures.append(h5['surfaceTemperatures']) 
        surfaceFractions.append(h5['surfaceFractions'] ) 
        LAI.append(h5['LAI'][()])
        windSpeed.append(5.0)
        windDirection.append(h5['windDirection10m'][()])
        n_absorbers.append(h5['n_absorbers'][()])
        landType.append(h5['landType'][()])
        soilType.append(h5['soilType'][()] ) 
        vegType.append(h5['vegType'][()] ) 
        waterType.append(h5['waterType'][()] )
        snowType.append(h5['snowType'][()] ) 
        iceType.append(h5['iceType'][()])
        storedTb.append( h5['Tb_atms'][0:22] )
        storedEmis.append( h5['emissivity_atms'][0:22] )
    print("Running Forward with 10 threads")
    start = time.time()
    forwardTb, forwardTransmission,\
    forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        np.asarray(zenithAngle).T, np.asarray(scanAngle).T,np.asarray(azimuthAngle).T, np.asarray(solarAngle).T, np.asarray(year), np.asarray(month), np.asarray(day), nChan, \
                        np.asarray(pressureLevels).T, np.asarray(pressureLayers).T, np.asarray(temperatureLayers).T, np.asarray(humidityLayers).T, np.asarray(ozoneConcLayers).T,\
                        np.asarray(co2ConcLayers).T,\
                        np.asarray(aerosolEffectiveRadius).T, np.asarray(aerosolConcentration).T, np.asarray(aerosolType).T, \
                        np.asarray(cloudEffectiveRadius).T, np.asarray(cloudConcentration).T, np.asarray(cloudType).T, np.asarray(cloudFraction).T, np.asarray(climatology).T, \
                        np.asarray(surfaceTemperatures).T, np.asarray(surfaceFractions).T, np.asarray(LAI), salinity*np.ones(len(LAI)), np.asarray(windSpeed).T, np.asarray(windDirection).T, np.asarray(n_absorbers).T,\
                        np.asarray(landType), np.asarray(soilType), np.asarray(vegType), np.asarray(waterType), np.asarray(snowType), np.asarray(iceType), 10)

    print("Running K Matrix with 10 threads.")
    kTb, kTransmission, temperatureJacobian, humidityJacobian, ozoneJacobian,\
    kEmissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        np.asarray(zenithAngle).T, np.asarray(scanAngle).T,np.asarray(azimuthAngle).T, np.asarray(solarAngle).T,np.asarray( year), np.asarray(month), np.asarray(day),nChan, \
                        np.asarray(pressureLevels).T, np.asarray(pressureLayers).T, np.asarray(temperatureLayers).T, np.asarray(humidityLayers).T, np.asarray(ozoneConcLayers).T,\
                        np.asarray(co2ConcLayers).T,\
                        np.asarray(aerosolEffectiveRadius).T, np.asarray(aerosolConcentration).T, np.asarray(aerosolType).T, \
                        np.asarray(cloudEffectiveRadius).T, np.asarray(cloudConcentration).T, np.asarray(cloudType).T, np.asarray(cloudFraction).T, np.asarray(climatology).T, \
                        np.asarray(surfaceTemperatures).T, np.asarray(surfaceFractions).T, np.asarray(LAI), salinity*np.ones(len(LAI)), np.asarray(windSpeed).T, np.asarray(windDirection).T, np.asarray(n_absorbers).T,\
                        np.asarray(landType), np.asarray(soilType), np.asarray(vegType), np.asarray(waterType), np.asarray(snowType), np.asarray(iceType), 10)
 

    end = time.time()
    print('pycrtm took',end-start)

    if ( all( np.abs( forwardTb.flatten() - np.asarray(storedTb).T.flatten() ) <= 1e-5)  and all( np.abs( kTb.flatten() - np.asarray(storedTb).T.flatten() ) <= 1e-5) ):
        print("Yay! all values are close enough to what CRTM test program produced!")
    else: 
        print("Boo! something failed. Look at all_200 plots")
        wavenumbers = np.linspace(1,23,22)
        plt.figure()
        plt.plot(wavenumbers,forwardTb-np.asarray(storedTb).T ) 
        plt.savefig(os.path.join(thisDir,'all_200'+'_spectrum_forward.png'))
        plt.figure()
        plt.plot(wavenumbers,forwardEmissivity-np.asarray(storedEmis).T)
        plt.savefig(os.path.join(thisDir,'all_200'+'_emissivity_forward.png')) 
    
        wavenumbers = np.linspace(1,23,22)
        plt.figure()
        plt.plot(wavenumbers,kTb-np.asarray(storedTb).T)
        plt.savefig(os.path.join(thisDir,'all_200'+'_spectrum_k.png'))
        plt.figure()
        plt.plot(wavenumbers,kEmissivity-np.asarray(storedEmis).T)
        plt.savefig(os.path.join(thisDir,'all_200'+'_emissivity_k.png')) 
        sys.exit("Boo! didn't pass tolerance with CRTM test program.")


if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(parentDir,'crtm.cfg') ) 
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'atms_npp'
    main(coefficientPath, sensor_id)
 
