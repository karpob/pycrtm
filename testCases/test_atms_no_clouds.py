#!/usr/bin/env python3
import configparser
import os, h5py, sys 
import numpy as np
from matplotlib import pyplot as plt
thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(thisDir)
sys.path.insert(0,parentDir)
from pycrtm import pycrtm
 
def main(coefficientPath, sensor_id):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    salinity = 35.0
    for c in cases:
        h5 = h5py.File(os.path.join(thisDir, 'data',c) , 'r')
        nChan = 22
        nans = -9999.9* np.ones(np.asarray(h5['aerosolConcentration']).shape)
        forwardTb, forwardTransmission,\
        forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, -10000.0, h5['solarAngle'].value, nChan, \
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], nans, h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], nans, h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, salinity, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)

        kTb, kTransmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, kEmissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, -10000.0, h5['solarAngle'].value, nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], nans, h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], nans, h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, salinity, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)
        

        wavenumbers = np.arange(22)
        diffK = kTb-h5['Tb_atms'][0:22]
        diffKemis = kEmissivity-h5['emissivity_atms'][0:22]
        
        if ( all(np.abs(diffKemis) <= 1e-10)  and all(np.abs(diffK) <= 1e-10) ):
            print ("Yay! we duplicated results from CRTM test program!")
        else:
            plt.figure()
            plt.plot(wavenumbers,kTb-h5['Tb_atms'][0:22])
            plt.savefig(os.path.join(thisDir,c+'_spectrum_k_matrix.png'))

            plt.figure()
            plt.plot(wavenumbers,kEmissivity-h5['emissivity_atms'][0:22])
            plt.savefig(os.path.join(thisDir,c+'_emissivity_k_matrix.png')) 

            plt.figure()
            plt.plot(wavenumbers,forwardTb-h5['Tb_atms'][0:22])
            plt.savefig(os.path.join(thisDir,c+'_spectrum_forward.png'))

            plt.figure()
            plt.plot(wavenumbers,forwardEmissivity-h5['emissivity_atms'][0:22])
            plt.savefig(os.path.join(thisDir,c+'_emissivity_forward.png')) 
            print("Boo! {} failed to pass a test. look at plots in {} for details.".format(c,thisDir) )

if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(parentDir,'crtm.cfg') ) 
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'atms_npp'
    main(coefficientPath, sensor_id)
 
