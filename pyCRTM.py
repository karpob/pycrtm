#!/usr/bin/env python3
import configparser
import os, h5py, sys 
import numpy as np
from matplotlib import pyplot as plt
thisDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,thisDir)
from pycrtm import pycrtm
from crtm_io import readSpcCoeff
from collections import namedtuple

def profilesCreate( nProfiles, nLevels, nAerosols=1, nClouds=1 ):
    keys  = [ 'P', 'T', 'Q', 'CO2', 'O3'] 
    p = {}
    for k in list(keys):
        p[k] = np.nan*np.ones([nProfiles,nLevels])

    p['Pi'] = np.nan*np.ones([nProfiles, nLevels+1])
    # satzen, sataz, sunzen, sunaz, scanangle
    p['Angles'] = np.nan*np.ones([nProfiles, 5])
    # P (2meter), T (2meter), Q (2meter), U(10meter), V(10meter), fetch 
    p['S2m'] = np.nan*np.ones([nProfiles,6])
    # skin T, salinity, snow_fraction, foam_fraction, fastem coef 1:5
    p['Skin'] = np.nan*np.zeros([nProfiles,10])
    # surftype, water type
    p['SurfType'] = np.nan*np.zeros([nProfiles,2])
    # latitude, longitude, elevation 
    p['SurfGeom'] = np.nan*np.zeros([nProfiles,3])
    # yy, mm, dd, hh, mm, ss
    p['DateTimes'] = np.zeros([nProfiles,6], dtype=int) 
    p['DateTimes'][:,0] = 2001
    p['DateTimes'][:,1] = 1
    p['DateTimes'][:,2] = 1
    # concentration, effective radius
    p['aerosols'] = np.nan*np.ones([nProfiles, nLevels, nAerosols, 2])
    p['aerosolType'] =-1 *np.ones([nProfiles,nAerosols], dtype =int)
    
    # concentration, effective radius
    p['clouds'] =  np.nan*np.ones([nProfiles, nLevels,  nClouds, 2])
    p['cloudType'] = -1 *np.ones([nProfiles,nClouds], dtype =int)
    p['cloudFraction'] = np.zeros([nProfiles,nLevels])

    p['LAI'] = np.zeros([nProfiles])
    # set so clouds are off unless you add something
    p['clouds'][:,:,0] = -1
    # surface 
    p['surfaceTemperatures'] = np.zeros([nProfiles,4])
    p['surfaceFractions'] = np.zeros([nProfiles,4])
    # land, soil, veg, water, snow, ice
    p['surfaceTypes'] = np.zeros([nProfiles,6], dtype=int)
    p['climatology'] = 6*np.ones([nProfiles], dtype=int)  # use usstd as default climatology for unspecified layers to 0.005 mbar (crtm will fill in the gaps if the user doesn't)

    p['windSpeed10m'] = np.zeros([nProfiles])
    p['windDirection10m'] = np.zeros([nProfiles])
    p['n_absorbers'] = np.zeros([nProfiles])

    profiles = namedtuple("Profiles", p.keys())(*p.values())
    return profiles

class pyCRTM:
    def __init__(self):
        self.coefficientPath = ''
        self.sensor_id = ''
        self.profiles = []
        self.Bt = []
        self.TauLevels = []
        self.surfEmisRefl = []
        self.TK = []
        self.QK = []
        self.O3K = []
        self.Wavenumbers = []
        self.wavenumbers = []
        self.wavenumber = []
        self.Wavenumber = []
        self.frequencyGHz = []
        self.wavelengthMicrons = []
        self.wavelengthCm = []
        self.nChan = 0
        self.nThreads = 1
    def loadInst(self):
        if ( os.path.exists( os.path.join(self.coefficientPath, self.sensor_id+'.SpcCoeff.bin') ) ):
            o = readSpcCoeff(os.path.join(self.coefficientPath, self.sensor_id+'.SpcCoeff.bin'))
            self.nChan = o['n_Channels']
            # For those who care to associate channel number with something physical:
            # just to save sanity put the permutations of (W/w)avenumber(/s) in here so things just go.
            self.wavenumbers = np.asarray(o['Wavenumber'])
            self.wavenumber = self.wavenumbers
            self.Wavenumber = self.wavenumbers 
            self.Wavenumbers = self.wavenumbers
            #For those more microwave oriented:
            self.frequencyGHz = 29.9792458 * self.wavenumbers
            self.wavelengthCm = 1.0/self.wavenumbers
            # And those who aren't interferometer oriented (people who like um): 
            self.wavelengthMicrons = 10000.0/self.wavenumbers
        else:
            print("Warning! {} doesn't exist!".format( os.path.join(self.coefficientPath, self.sensor_id+'.SpcCoeff.bin') ) )        
    def runDirect(self):
        items =dir(self.profiles) 
        #print(pycrtm.__doc__) 
        self.Bt, layerOpticalDepths,\
        self.surfEmisRefl  = pycrtm.wrap_forward( self.coefficientPath, self.sensor_id,\
                        self.profiles.Angles[:,0], self.profiles.Angles[:,4], self.profiles.Angles[:,1], self.profiles.Angles[:,2:4], self.profiles.DateTimes[:,0], self.profiles.DateTimes[:,1],self.profiles.DateTimes[:,2], self.nChan, \
                        self.profiles.Pi, self.profiles.P, self.profiles.T, self.profiles.Q, self.profiles.O3,\
                        self.profiles.CO2,\
                        self.profiles.aerosols[:,:,:,1], self.profiles.aerosols[:,:,:,0], self.profiles.aerosolType, \
                        self.profiles.clouds[:,:,:,1], self.profiles.clouds[:,:,:,0], self.profiles.cloudType, self.profiles.cloudFraction, self.profiles.climatology, \
                        self.profiles.surfaceTemperatures, self.profiles.surfaceFractions, self.profiles.LAI, self.profiles.S2m[:,1], self.profiles.windSpeed10m, self.profiles.windDirection10m, self.profiles.n_absorbers,\
                        self.profiles.surfaceTypes[:,0], self.profiles.surfaceTypes[:,1], self.profiles.surfaceTypes[:,2], self.profiles.surfaceTypes[:,3], self.profiles.surfaceTypes[:,4], self.profiles.surfaceTypes[:,5], self.nThreads )
        self.TauLevels = np.zeros(layerOpticalDepths.shape)
        nprofile, nchan, nlay = layerOpticalDepths.shape
        # should use python threading here!
        for p in list(range(nprofile)):
            for c in list(range(nchan)):
                self.TauLevels[p,c,:] = np.exp(-1.0*np.cumsum(layerOpticalDepths[p,c,:] ))
    def runK(self):
        self.Bt, layerOpticalDepths, self.TK, self.QK, self.O3K,\
        self.surfEmisRefl =  pycrtm.wrap_k_matrix(  self.coefficientPath, self.sensor_id,\
                        self.profiles.Angles[:,0], self.profiles.Angles[:,4], self.profiles.Angles[:,1], self.profiles.Angles[:,2:4], self.profiles.DateTimes[:,0], self.profiles.DateTimes[:,1],self.profiles.DateTimes[:,2], self.nChan, \
                        self.profiles.Pi, self.profiles.P, self.profiles.T, self.profiles.Q, self.profiles.O3,\
                        self.profiles.CO2,\
                        self.profiles.aerosols[:,:,:,1], self.profiles.aerosols[:,:,:,0], self.profiles.aerosolType, \
                        self.profiles.clouds[:,:,:,1], self.profiles.clouds[:,:,:,0], self.profiles.cloudType, self.profiles.cloudFraction, self.profiles.climatology, \
                        self.profiles.surfaceTemperatures, self.profiles.surfaceFractions, self.profiles.LAI, self.profiles.S2m[:,1], self.profiles.windSpeed10m, self.profiles.windDirection10m, self.profiles.n_absorbers,\
                        self.profiles.surfaceTypes[:,0], self.profiles.surfaceTypes[:,1], self.profiles.surfaceTypes[:,2], self.profiles.surfaceTypes[:,3], self.profiles.surfaceTypes[:,4], self.profiles.surfaceTypes[:,5], self.nThreads )

        self.TauLevels = np.zeros(layerOpticalDepths.shape)
        nprofile, nchan, nlay = layerOpticalDepths.shape
        # should use python threading here!
        for p in list(range(nprofile)):
            for c in list(range(nchan)):
                self.TauLevels[p,c,:] = np.exp(-1.0*np.cumsum(layerOpticalDepths[p,c,:]))


if __name__ == "__main__":
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'testCases','data') )
    cases.sort()
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(thisDir,'crtm.cfg') ) 


    profiles = profilesCreate( 4, 92 )
    
    for i,c in enumerate(cases):
        h5 = h5py.File(os.path.join(thisDir, 'testCases','data',c) , 'r')
        profiles.Angles[i,0] = h5['zenithAngle'][()]
        profiles.Angles[i,1] = 999.9 
        profiles.Angles[i,2] = 100.0  # 100 degrees zenith below horizon.
        profiles.Angles[i,3] = 0.0 # zero solar azimuth 
        profiles.Angles[i,4] = h5['scanAngle'][()]
        profiles.DateTimes[i,0] = 2001
        profiles.DateTimes[i,1] = 1
        profiles.DateTimes[i,2] = 1
        profiles.Pi[i,:] = np.asarray(h5['pressureLevels'] )
        profiles.P[i,:] = np.asarray(h5['pressureLayers'][()])
        profiles.T[i,:] = np.asarray(h5['temperatureLayers'])
        profiles.Q[i,:] = np.asarray(h5['humidityLayers'])
        profiles.O3[i,:] = np.asarray(h5['ozoneConcLayers'])
        profiles.CO2[i,:] = np.asarray(h5['co2ConcLayers'])
        profiles.clouds[i,:,0,0] = np.asarray(h5['cloudConcentration'])
        profiles.clouds[i,:,0,1] = np.asarray(h5['cloudEffectiveRadius'])
        profiles.aerosols[i,:,0,0] = np.asarray(h5['aerosolConcentration'])
        profiles.aerosols[i,:,0,1] = np.asarray(h5['aerosolEffectiveRadius'])
        profiles.aerosolType[i] = h5['aerosolType'][()]
        profiles.cloudType[i] = h5['cloudType'][()]
        profiles.cloudFraction[i,:] = h5['cloudFraction'][()]
        profiles.climatology[i] = h5['climatology'][()]
        profiles.surfaceFractions[i,:] = h5['surfaceFractions']
        profiles.surfaceTemperatures[i,:] = h5['surfaceTemperatures']
        profiles.S2m[i,1] = 33.0 # just use salinity out of S2m for the moment.
        profiles.windSpeed10m[i] = 5.0
        profiles.windDirection10m[i] = h5['windDirection10m'][()]
        profiles.n_absorbers[i] = h5['n_absorbers'][()]
        # land, soil, veg, water, snow, ice
        profiles.surfaceTypes[i,0] = h5['landType'][()]
        profiles.surfaceTypes[i,1] = h5['soilType'][()]
        profiles.surfaceTypes[i,2] = h5['vegType'][()]
        profiles.surfaceTypes[i,3] = h5['waterType'][()]
        profiles.surfaceTypes[i,4] = h5['snowType'][()]
        profiles.surfaceTypes[i,5] = h5['iceType'][()]
        profiles.LAI[i] = h5['LAI'][()]
        h5.close()
    crtmOb = pyCRTM()
    crtmOb.profiles = profiles
    crtmOb.coefficientPath = pathInfo['CRTM']['coeffs_dir']
    crtmOb.sensor_id = 'atms_npp'
    crtmOb.nThreads = 4
    crtmOb.loadInst()
    crtmOb.runDirect()
    crtmOb.runK()




