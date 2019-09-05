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
# Absorber IDs taken from CRTM.
gases = {}
gases['Q']     = 1  # H2O for anyone not NWP focused ;)
gases['CO2']   = 2
gases['O3']    = 3
gases['N2O']   = 4
gases['CO']    = 5
gases['CH4']   = 6
gases['O2']    = 7
gases['NO']    = 8
gases['SO2']   = 9
gases['NO2']   = 10
gases['NH3']   = 11
gases['HNO3']  = 12
gases['OH']    = 13
gases['HF']    = 14
gases['HCl']   = 15
gases['HBr']   = 16
gases['HI']    = 17
gases['ClO']   = 18
gases['OCS']   = 19
gases['H2CO']  = 20
gases['HOCl']  = 21
gases['N2']    = 22
gases['HCN']   = 23
gases['CH3l']  = 24
gases['H2O2']  = 25
gases['C2H2']  = 26
gases['C2H6']  = 27
gases['PH3']   = 28
gases['COF2']  = 29
gases['SF6']   = 30
gases['H2S']   = 31
gases['HCOOH'] = 32
def profilesCreate( nProfiles, nLevels, nAerosols=1, nClouds=1, additionalGases=[] ):
    keys  = [ 'P', 'T', 'Q', 'O3']
    for g in additionalGases:
        if g in list(gases.keys()) and g not in keys:
            keys.append(g)
        elif g == 'H2O' or g.lower() == 'water' or g=='ozone':
            print("You worry too much, of course we have {}! Water and Ozone are always turned on.".format(g))
        else: 
            print("Warning! I don't know this gas: {}! I can't add it to the simulation!".format(g))
            print("You could pick one of these instead:")
            for gg in list(gases.keys()): print(gg)

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

def apply_avg(Pout, Pin, Xin):
    """
    Use LayerAvg function from CRTM to Interpolate User pressure Layers to Coefficient levels. 
    """
    Xout = pycrtm.applyavg(np.log(Pout),np.log(Pin),Xin)
    return Xout

class pyCRTM:
    def __init__(self):
        self.coefficientPath = ''
        self.sensor_id = ''
        self.profiles = []
        self.traceConc = []
        self.traceIds = []
        self.usedGases = []
        self.Bt = []
        self.TauLevels = []
        self.surfEmisRefl = []
        self.TK = []
        self.QK = []
        self.O3K = []
        self.CO2K = []
        self.N2OK = []
        self.CH4K = []
        self.COK = []
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
            self.wmo_sensor_id = o['wmo_sensor_id']
            self.wmo_satellite_id = o['wmo_satellite_id']
        else:
            print("Warning! {} doesn't exist!".format( os.path.join(self.coefficientPath, self.sensor_id+'.SpcCoeff.bin') ) )        
    def setupGases(self):
        max_abs = int(max(self.profiles.n_absorbers))
        nprof, nlay = self.profiles.T.shape 
        self.traceConc = np.zeros([nprof,nlay,max_abs])
        self.traceIds = np.zeros(max_abs, dtype=np.int)
        availableGases = list(gases.keys())
        profileItems = list(self.profiles._asdict().keys())
       
        for p in profileItems:
            if (p in availableGases):self.usedGases.append(p)
        
        for i,g in enumerate(self.usedGases):
            self.traceConc[:,:,i] = self.profiles._asdict()[g][:,:]
            self.traceIds[i] = gases[g]
        
    def runDirect(self):
        
        self.setupGases() 
        items =dir(self.profiles) 
        #print(pycrtm.__doc__) 
        self.Bt, layerOpticalDepths,\
        self.surfEmisRefl  = pycrtm.wrap_forward( self.coefficientPath, self.sensor_id,\
                        self.profiles.Angles[:,0], self.profiles.Angles[:,4], self.profiles.Angles[:,1], self.profiles.Angles[:,2:4], self.profiles.DateTimes[:,0], self.profiles.DateTimes[:,1],self.profiles.DateTimes[:,2], self.nChan, \
                        self.profiles.Pi, self.profiles.P, self.profiles.T, self.traceConc,self.traceIds,\
                        self.profiles.aerosols[:,:,:,1], self.profiles.aerosols[:,:,:,0], self.profiles.aerosolType, \
                        self.profiles.clouds[:,:,:,1], self.profiles.clouds[:,:,:,0], self.profiles.cloudType, self.profiles.cloudFraction, self.profiles.climatology, \
                        self.profiles.surfaceTemperatures, self.profiles.surfaceFractions, self.profiles.LAI, self.profiles.S2m[:,1], self.profiles.windSpeed10m, self.profiles.windDirection10m, self.profiles.n_absorbers,\
                        self.profiles.surfaceTypes[:,0], self.profiles.surfaceTypes[:,1], self.profiles.surfaceTypes[:,2], self.profiles.surfaceTypes[:,3], self.profiles.surfaceTypes[:,4], self.profiles.surfaceTypes[:,5], self.nThreads )
        self.TauLevels = np.zeros(layerOpticalDepths.shape)
        nprofile, nchan, nlay = layerOpticalDepths.shape
        # should use python threading here!
        # TauLevels following RTTOV convention: 
        # "Transmittance from each user pressure level to TOA." 
        for p in list(range(nprofile)):
            for c in list(range(nchan)):
                self.TauLevels[p,c,:] = np.exp(-1.0*np.cumsum(layerOpticalDepths[p,c,:] ))
    def runK(self):
        self.setupGases() 
             
        self.Bt, layerOpticalDepths, self.TK, traceK,\
        self.surfEmisRefl =  pycrtm.wrap_k_matrix(  self.coefficientPath, self.sensor_id,\
                        self.profiles.Angles[:,0], self.profiles.Angles[:,4], self.profiles.Angles[:,1], self.profiles.Angles[:,2:4], self.profiles.DateTimes[:,0], self.profiles.DateTimes[:,1],self.profiles.DateTimes[:,2], self.nChan, \
                        self.profiles.Pi, self.profiles.P, self.profiles.T, \
                        self.traceConc, self.traceIds,\
                        self.profiles.aerosols[:,:,:,1], self.profiles.aerosols[:,:,:,0], self.profiles.aerosolType, \
                        self.profiles.clouds[:,:,:,1], self.profiles.clouds[:,:,:,0], self.profiles.cloudType, self.profiles.cloudFraction, self.profiles.climatology, \
                        self.profiles.surfaceTemperatures, self.profiles.surfaceFractions, self.profiles.LAI, self.profiles.S2m[:,1], self.profiles.windSpeed10m, self.profiles.windDirection10m, self.profiles.n_absorbers,\
                        self.profiles.surfaceTypes[:,0], self.profiles.surfaceTypes[:,1], self.profiles.surfaceTypes[:,2], self.profiles.surfaceTypes[:,3], self.profiles.surfaceTypes[:,4], self.profiles.surfaceTypes[:,5], self.nThreads )
        for i,ids in enumerate(list(self.traceIds)):
            # I think I can do something smarter here in python to contruct self.QK etc through an execute, or something along those lines?
            if(ids == gases['Q']):   self.QK   = traceK[:,:,:,i]
            if(ids == gases['O3']):  self.O3K  = traceK[:,:,:,i]
            if(ids == gases['CH4']): self.CH4K = traceK[:,:,:,i]
            if(ids == gases['CO2']): self.CO2K = traceK[:,:,:,i]
            if(ids == gases['CO']):  self.COK  = traceK[:,:,:,i]
            if(ids == gases['N2O']): self.N2OK = traceK[:,:,:,i]
        # if we don't have any "weird" gases, empty out traceK,traceConc to save on RAM.
        if not any(g in self.usedGases for g in  ['Q', 'O3', 'CH4', 'CO','CO2', 'N2O']):
            print("saving on RAM")
            self.traceK = []
            self.traceConc = []     
      

        self.TauLevels = np.zeros(layerOpticalDepths.shape)
        nprofile, nchan, nlay = layerOpticalDepths.shape
        # should use python threading here!
        # TauLevels following RTTOV convention: 
        # "Transmittance from each user pressure level to TOA." 
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




