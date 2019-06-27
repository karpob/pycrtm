#!/usr/bin/env python3
import configparser
import os, h5py 
from pycrtm import pycrtm 
import numpy as np
from matplotlib import pyplot as plt

def main(coefficientPath, sensor_id):
    cases = os.listdir('testCases/')
    salinity = 35.0
    for c in cases:
        h5 = h5py.File(os.path.join('testCases',c) , 'r')
        nChan = 22
        forwardTb, forwardTransmission,\
        forwardEmissivity = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, h5['azimuthAngle'].value, h5['solarAngle'].value, nChan, \
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
                        h5['zenithAngle'].value, h5['scanAngle'].value, h5['azimuthAngle'].value, h5['solarAngle'].value, nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, salinity, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)
        

        wavenumbers = np.arange(22)
 
        plt.figure()
        plt.plot(wavenumbers,kTb-h5['Tb_atms'][0:22])
        plt.savefig(c+'_spectrum_k_matrix.png')

        plt.figure()
        plt.plot(wavenumbers,kEmissivity-h5['emissivity_atms'][0:22])
        plt.savefig(c+'_emissivity_k_matrix.png') 

        plt.figure()
        plt.plot(wavenumbers,forwardTb-h5['Tb_atms'][0:22])
        plt.savefig(c+'_spectrum_forward.png')

        plt.figure()
        plt.plot(wavenumbers,forwardEmissivity-h5['emissivity_atms'][0:22])
        plt.savefig(c+'_emissivity_forward.png') 
    print('done.')
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read('crtm.cfg')
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'atms_npp'
    #sensor_id = 'cris_npp' 
    main(coefficientPath, sensor_id)
 
