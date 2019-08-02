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
        nans = -9999.9* np.ones(np.asarray(h5['aerosolConcentration']).shape)
        forwardTb, forwardTransmission,\
        forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'][()], h5['scanAngle'][()], 999.9, h5['solarAngle'][()], nChan, \
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], nans, h5['aerosolType'][()], \
                        h5['cloudEffectiveRadius'], nans, h5['cloudType'][()], h5['cloudFraction'], h5['climatology'][()], \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'][()], salinity, h5['windSpeed10m'][()], h5['windDirection10m'][()], h5['n_absorbers'][()],\
                        h5['landType'][()], h5['soilType'][()], h5['vegType'][()], h5['waterType'][()], h5['snowType'][()], h5['iceType'][()], 1)

        kTb, kTransmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, kEmissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        h5['zenithAngle'][()], h5['scanAngle'][()], 999.9, h5['solarAngle'][()], nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], nans, h5['aerosolType'][()], \
                        h5['cloudEffectiveRadius'], nans, h5['cloudType'][()], h5['cloudFraction'], h5['climatology'][()], \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'][()], salinity, h5['windSpeed10m'][()], h5['windDirection10m'][()], h5['n_absorbers'][()],\
                        h5['landType'][()], h5['soilType'][()], h5['vegType'][()], h5['waterType'][()], h5['snowType'][()], h5['iceType'][()],1)
        

        
        diffK = kTb.flatten()-h5['Tb']
        diffKemis = kEmissivity.flatten()-h5['emissivity']
        
        if ( all(np.abs(diffKemis) <= 1e-8)  and all(np.abs(diffK) <= 1e-8) ):
            print ("Yay! we duplicated results from CRTM test program!")
        else:
            h5wav = h5py.File(os.path.join(thisDir,'cris_wavenumbers.h5'),'r')
            wavenumbers = np.asarray(h5wav['wavenumbers'])
            
            plt.figure()
            plt.plot(wavenumbers,diffK)
            plt.savefig( os.path.join(thisDir,c+'_spectrum_k_matrix.png') )

            plt.figure()
            plt.plot(wavenumbers,diffKemis)
            plt.savefig( os.path.join(thisDir,c+'_emissivity_k_matrix.png') ) 

            plt.figure()
            plt.plot(wavenumbers,forwardTb.flatten()-h5['Tb'])
            plt.savefig( os.path.join(thisDir,c+'_spectrum_forward.png') )

            plt.figure()
            plt.plot(wavenumbers,forwardEmissivity.flatten()-h5['emissivity'])
            plt.savefig( os.path.join(thisDir,c+'_emissivity_forward.png') ) 
            print("Boo! {} failed to pass a test. look at plots for details in {}.".format(c,thisDir))

if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join( parentDir, 'crtm.cfg' ) )
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'cris_npp' 
    main(coefficientPath, sensor_id)
 
