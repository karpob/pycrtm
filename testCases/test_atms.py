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
        forwardTb, forwardTransmission,\
        forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'][()], h5['scanAngle'][()], 999.9, h5['solarAngle'][()], nChan, \
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'][()], \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'][()], h5['cloudFraction'], h5['climatology'][()], \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'][()], salinity, h5['windSpeed10m'][()], h5['windDirection10m'][()], h5['n_absorbers'][()],\
                        h5['landType'][()], h5['soilType'][()], h5['vegType'][()], h5['waterType'][()], h5['snowType'][()], h5['iceType'][()])

        kTb, kTransmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, kEmissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        h5['zenithAngle'][()], h5['scanAngle'][()], 999.9, h5['solarAngle'][()], nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'][()], \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'][()], h5['cloudFraction'], h5['climatology'][()], \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'][()], salinity, h5['windSpeed10m'][()], h5['windDirection10m'][()], h5['n_absorbers'][()],\
                        h5['landType'][()], h5['soilType'][()], h5['vegType'][()], h5['waterType'][()], h5['snowType'][()], h5['iceType'][()])
        

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
            sys.exit("Boo! {} failed to pass a test. look at plots in {} for details.".format(c,thisDir) )

if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(parentDir,'crtm.cfg') ) 
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'atms_npp'
    main(coefficientPath, sensor_id)
 
