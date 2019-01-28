#!/usr/bin/env python3
import configparser
from pycrtm import pycrtm 
import numpy as np
def main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m):
    p = {}
    p['US_Std'] ={}
    p['Trop'] = {} 

    p['US_Std']['pressureLevels'],\
    p['US_Std']['pressureLayers'],\
    p['US_Std']['temperatureLayers'],\
    p['US_Std']['humidityLayers'],\
    p['US_Std']['ozoneConcLayers'] = pycrtm.test_data_us_std()

    p['Trop']['pressureLevels'],\
    p['Trop']['pressureLayers'],\
    p['Trop']['temperatureLayers'],\
    p['Trop']['humidityLayers'],\
    p['Trop']['ozoneConcLayers'] = pycrtm.test_data_tropical()

    for k in list(p.keys()):
        N_LAYERS = p[k]['pressureLayers'].shape[0]
        
        Tb = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                              zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan, \
                              p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                              surfaceType, surfaceTemperature, windSpeed10m, windDirection10m )
        print(Tb)
        
        Tb, Transmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan,\
                        p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                        surfaceType, surfaceTemperature, windSpeed10m, windDirection10m )
        print(ozoneJacobian)
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read('crtm.cfg')
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    #sensor_id = 'atms_n20'
    sensor_id = 'cris_npp'
    zenithAngle = float(30.0)
    scanAngle = float(26.37293341421)
    azimuthAngle = float(0.0)
    solarAngle = float(0.0)
    nChan = int(1305)
    surfaceType = int(1)
    surfaceTemperature = float(293.15)
    windSpeed10m = float(0.0)
    windDirection10m = float(0.0)
    # assuming this is zenith? probably need Az too?
    main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m)
 
