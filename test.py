#!/usr/bin/env python3
import configparser
import h5py 
from pycrtm import pycrtm 
import numpy as np
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as mplcm
from cycler import cycler 
from matplotlib import pyplot as plt

def plotJacobians(chan_list, p, t, q, jacobians, instrument, wavenumbers, ofWhat):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    jacSub = []
    if(ofWhat.lower()=='ozone'):
        for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*q
    elif(ofWhat.lower()=='water'):
         for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*q
    elif(ofWhat.lower()=='temperature'):
          for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*t

    for i in range(jacobians.shape[0]):
        if(i+1 in chan_list):
            plt.plot( jacobians[i,:], p )
            jacSub.append(jacobians[i,:])
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wavenumbers[c-1], 10000.0/wavenumbers[c-1]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    if(ofWhat.lower() =='temperature'): plt.xlabel('Jacobian [K]')
    else: plt.xlabel('Jacobian [K]')
    plt.savefig(instrument+'_jacobians_'+ofWhat+'.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
    plt.close()

    plt.figure()
    if(ofWhat.lower() == 'temperature'): plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel')
    else: plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel')
    jacSub = np.asarray(jacSub)
    #plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = 0.3 ))#vmax = wf[:,0:99].max() ) )
    plt.pcolor(np.arange(len(chan_list)+1), p, jacSub[:,:].T ) 
    plt.colorbar()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Instrument Channel')
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xticks(np.arange(len(chan_list)), chan_list, rotation='vertical')
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    plt.tight_layout() 
    plt.savefig(instrument+'_jacobians_'+ofWhat+'_pcolor.png')
    plt.close()

def main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m):
    """    
    p = {}
    p['US_Std'] ={}
    with open('Temperature_CRTM1.bin') as f:
        p['US_Std']['temperatureLayers'] = np.fromfile(f,dtype='<f8')
    with open('aerosolConcentration_CRTM1.bin') as f:
        p['US_Std']['aerosolConcentration'] = np.fromfile(f,dtype='<f8')
    with open('aerosolEffectiveRadius_CRTM1.bin') as f:
        p['US_Std']['aerosolEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('cloudEffectiveRadius_CRTM1.bin') as f:
        p['US_Std']['cloudEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('co2_CRTM1.bin') as f:
        p['US_Std']['co2ConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('ozone_CRTM1.bin') as f:
        p['US_Std']['ozoneConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('pressureLevels_CRTM1.bin') as f:
        p['US_Std']['pressureLevels'] = np.fromfile(f,dtype='<f8')
    with open('pressure_CRTM1.bin') as f:
        p['US_Std']['pressureLayers'] = np.fromfile(f,dtype='<f8')
    with open('waterContent_CRTM1.bin') as f:
        p['US_Std']['cloudConcentration'] = np.fromfile(f,dtype='<f8')
    with open('waterVapor_CRTM1.bin') as f:
        p['US_Std']['humidityLayers'] = np.fromfile(f,dtype='<f8')
    with open('Temperature_CRTM2.bin') as f:
        p['US_Std']['temperatureLayers'] = np.fromfile(f,dtype='<f8')
    with open('aerosolConcentration_CRTM2.bin') as f:
        p['US_Std']['aerosolConcentration'] = np.fromfile(f,dtype='<f8')
    with open('aerosolEffectiveRadius_CRTM2.bin') as f:
        p['US_Std']['aerosolEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('cloudEffectiveRadius_CRTM2.bin') as f:
        p['US_Std']['cloudEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('co2_CRTM2.bin') as f:
        p['US_Std']['co2ConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('ozone_CRTM2.bin') as f:
        p['US_Std']['ozoneConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('pressureLevels_CRTM2.bin') as f:
        p['US_Std']['pressureLevels'] = np.fromfile(f,dtype='<f8')
    with open('pressure_CRTM2.bin') as f:
        p['US_Std']['pressureLayers'] = np.fromfile(f,dtype='<f8')
    with open('waterContent_CRTM2.bin') as f:
        p['US_Std']['cloudConcentration'] = np.fromfile(f,dtype='<f8')
    with open('waterVapor_CRTM2.bin') as f:
        p['US_Std']['humidityLayers'] = np.fromfile(f,dtype='<f8')
    p['US_Std']['cloudFraction'] = np.zeros( p['US_Std']['cloudConcentration'].shape)
    #p['US_Std']['cloudFraction'][np.where(p['US_Std']['cloudConcentration']> 0.0)] = 0.1426
    p['US_Std']['cloudFraction'][np.where(p['US_Std']['cloudConcentration']> 0.0)] = 0.0
    
    p['US_Std']['aerosolType'] = 1 # 1 Dust., 2 Sea salt
    p['US_Std']['cloudType'] = 1 # 1 Water ( I think ) 3, rain
    # Land, water, snow, ice 
    surfaceTemperatures = np.asarray([272.0, 275.0, 265.0, 269.0])
    #surfaceTemperatures = np.asarray([318.0, 0.0, 0.0, 0.0])
    surfaceFractions = np.asarray([0.1, 0.5, 0.25, 0.15]) 
    #surfaceFractions = np.asarray([1.0, 0.0, 0.0, 0.0]) 
    n_absorbers = 2
    climatology = 6
    LAI = 0.17
    #LAI = 0.65
    TUNDRA_SURFACE_TYPE         = 10  # NPOESS Land surface type for IR/VIS Land SfcOptics
    SCRUB_SURFACE_TYPE          =  7  # NPOESS Land surface type for IR/VIS Land SfcOptics
    COARSE_SOIL_TYPE            =  1  # Soil type                for MW land SfcOptics
    GROUNDCOVER_VEGETATION_TYPE =  7  # Vegetation type          for MW Land SfcOptics
    BARE_SOIL_VEGETATION_TYPE   = 11  # Vegetation type          for MW Land SfcOptics
    SEA_WATER_TYPE              =  1  # Water type               for all SfcOptics
    FRESH_SNOW_TYPE             =  2  # NPOESS Snow type         for IR/VIS SfcOptics
    FRESH_ICE_TYPE              =  1  # NPOESS Ice type          for IR/VIS SfcOptics


    landType = TUNDRA_SURFACE_TYPE
    soilType = COARSE_SOIL_TYPE
    vegType = GROUNDCOVER_VEGETATION_TYPE
    waterType = SEA_WATER_TYPE
    snowType = FRESH_SNOW_TYPE
    iceType = FRESH_ICE_TYPE

    h5 = h5py.File('case1.h5','w')
    for k in list(p['US_Std'].keys()):
        h5.create_dataset(k,data = p['US_Std'][k])
    h5.create_dataset('landType',data = landType)
    h5.create_dataset('soilType',data = soilType)
    h5.create_dataset('vegType',data = vegType)
    h5.create_dataset('waterType',data = waterType)
    h5.create_dataset('snowType', data = snowType)
    h5.create_dataset('iceType', data = iceType)
    h5.create_dataset('LAI', data = LAI)
    h5.create_dataset('n_absorbers', data = 2)
    h5.create_dataset('climatology', data = climatology)
    h5.create_dataset('surfaceTemperatures', data = surfaceTemperatures)
    h5.create_dataset('surfaceFractions', data = surfaceFractions)
    h5.create_dataset('zenithAngle', data = zenithAngle)
    h5.create_dataset('scanAngle', data = scanAngle)
    h5.create_dataset('azimuthAngle', data = azimuthAngle)
    h5.create_dataset('solarAngle', data = solarAngle)
    h5.create_dataset('windDirection10m', data = windDirection10m)
    h5.create_dataset('windSpeed10m', data = windDirection10m)
    """
    h5 = h5py.File('case4.h5','r')

    chan_list = [577, 607, 626, 650, 667]
    for k in [1,]:
        print(list(h5.keys()))
        Tb1, Transmission1,\
        emissivity1 = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, h5['azimuthAngle'].value, h5['solarAngle'].value, nChan, \
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)

        print('done forward.')

        Tb, Transmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, emissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        h5['zenithAngle'].value, h5['scanAngle'].value, h5['azimuthAngle'].value, h5['solarAngle'].value, nChan,\
                        h5['pressureLevels'], h5['pressureLayers'], h5['temperatureLayers'], h5['humidityLayers'], h5['ozoneConcLayers'],\
                        h5['co2ConcLayers'],\
                        h5['aerosolEffectiveRadius'], h5['aerosolConcentration'], h5['aerosolType'].value, \
                        h5['cloudEffectiveRadius'], h5['cloudConcentration'], h5['cloudType'].value, h5['cloudFraction'], h5['climatology'].value, \
                        h5['surfaceTemperatures'], h5['surfaceFractions'], h5['LAI'].value, h5['windSpeed10m'].value, h5['windDirection10m'].value, h5['n_absorbers'].value,\
                        h5['landType'].value, h5['soilType'].value, h5['vegType'].value, h5['waterType'].value, h5['snowType'].value, h5['iceType'].value)

        print('done k matrix.')
        h5wav = h5py.File('cris_wavenumbers.h5','r')
        wavenumbers = np.asarray(h5wav['wavenumbers'])
        plotJacobians(chan_list, np.asarray(h5['pressureLayers']), np.asarray(h5['temperatureLayers']), np.asarray(h5['ozoneConcLayers']), ozoneJacobian, 'CrIS_'+str(k)+'_', wavenumbers, 'ozone')
        plotJacobians(chan_list, np.asarray(h5['pressureLayers']), np.asarray(h5['temperatureLayers']), np.asarray(h5['ozoneConcLayers']), temperatureJacobian, 'CrIS_'+str(k)+'_', wavenumbers, 'temperature')
        plotJacobians(chan_list, np.asarray(h5['pressureLayers']), np.asarray(h5['temperatureLayers']), np.asarray(h5['humidityLayers']), humidityJacobian, 'CrIS_'+str(k)+'_', wavenumbers, 'water')
        plt.figure()
        plt.plot(wavenumbers,Tb-h5['Tb'])
        plt.savefig('spectrum.png')

        plt.figure()
        plt.plot(wavenumbers,emissivity-h5['emissivity'])
        plt.savefig('emissivity.png') 
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
    windSpeed10m = float(0.0)
    windDirection10m = float(0.0)
    surfaceTemperature = 0
    # assuming this is zenith? probably need Az too?
    main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m)
 
