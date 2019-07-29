#!/usr/bin/env python3
import configparser
import os, h5py,sys 
import numpy as np
from matplotlib import pyplot as plt
thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(thisDir)
sys.path.insert(0,parentDir)
from pycrtm import pycrtm

def main(coefficientPath, sensor_id):
    salinity = 35.0
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') ) 
    cases.sort()
    for c in cases:
        h5 = h5py.File(os.path.join(thisDir, 'data', c) , 'r')
        nChan = np.asarray(h5['Tb']).shape[0] 
        forwardTb, forwardTransmission,\
        forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, -10000.0, h5['solarAngle'].value, nChan, \
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, salinity, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)

        kTb, kTransmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, kEmissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, -10000.0, h5['solarAngle'].value, nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, salinity, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)
        

        
        diffK = kTb-h5['Tb']
        diffKemis = kEmissivity-h5['emissivity']
        
        if ( all(np.abs(diffKemis) <= 0.0)  and all(np.abs(diffK) <= 0.0) ):
            print ("Yay! we duplicated results from CRTM test program!")
        else:
            h5wav = h5py.File(os.path.join(thisDir,'cris_wavenumbers.h5'),'r')
            wavenumbers = np.asarray(h5wav['wavenumbers'])
            
            plt.figure()
            plt.plot(wavenumbers,kTb-h5['Tb'])
            plt.savefig( os.path.join(thisDir,c+'_spectrum_k_matrix.png') )

            plt.figure()
            plt.plot(wavenumbers,kEmissivity-h5['emissivity'])
            plt.savefig( os.path.join(thisDir,c+'_emissivity_k_matrix.png') ) 

            plt.figure()
            plt.plot(wavenumbers,forwardTb-h5['Tb'])
            plt.savefig( os.path.join(thisDir,c+'_spectrum_forward.png') )

            plt.figure()
            plt.plot(wavenumbers,forwardEmissivity-h5['emissivity'])
            plt.savefig( os.path.join(thisDir,c+'_emissivity_forward.png') ) 
            sys.exit("Boo! {} failed to pass a test. look at plots for details in {}.".format(c,thisDir))

if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join( parentDir, 'crtm.cfg' ) )
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'cris_npp' 
    main(coefficientPath, sensor_id)
 
