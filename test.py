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
    p = {}
    p['US_Std'] ={}
    #p['Trop'] = {} 

    p['US_Std']['pressureLevels'],\
    p['US_Std']['pressureLayers'],\
    p['US_Std']['temperatureLayers'],\
    p['US_Std']['humidityLayers'],\
    p['US_Std']['ozoneConcLayers'],\
    p['US_Std']['aerosolEffectiveRadius'],\
    p['US_Std']['aerosolConcentration'],\
    p['US_Std']['aerosolType'] = pycrtm.test_data_us_std()

    #p['Trop']['pressureLevels'],\
    #p['Trop']['pressureLayers'],\
    #p['Trop']['temperatureLayers'],\
    #p['Trop']['humidityLayers'],\
    #p['Trop']['ozoneConcLayers'] = pycrtm.test_data_tropical()
    chan_list = [577, 607, 626, 650, 667]
    for k in list(p.keys()):
        N_LAYERS = p[k]['pressureLayers'].shape[0]
        """ 
        Tb,transmission = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                              zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan, \
                              p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                              p[k]['aerosolEffectiveRadius'], p[k]['aerosolConcentration'], p[k]['aerosolType'],\
                              surfaceType, surfaceTemperature, windSpeed10m, windDirection10m )
        print(Tb,transmission)
        """
        Tb, Transmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan,\
                        p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                        p[k]['aerosolEffectiveRadius'], p[k]['aerosolConcentration'], p[k]['aerosolType'],\
                        surfaceType, surfaceTemperature, windSpeed10m, windDirection10m )

        h5 = h5py.File('cris_wavenumbers.h5','r')
        wavenumbers = np.asarray(h5['wavenumbers'])
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['ozoneConcLayers'], ozoneJacobian, 'CrIS_'+k+'_', wavenumbers, 'ozone')
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['ozoneConcLayers'], temperatureJacobian, 'CrIS_'+k+'_', wavenumbers, 'temperature')
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], humidityJacobian, 'CrIS_'+k+'_', wavenumbers, 'water')
       
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read('crtm.cfg')
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    #sensor_id = 'atms_n20'
    sensor_id = 'cris_npp'
    zenithAngle = float(0) #float(30.0)
    scanAngle = float(0) #float(26.37293341421)
    azimuthAngle = float(0.0)
    solarAngle = float(0.0)
    nChan = int(1305)
    surfaceType = int(1)
    surfaceTemperature = float(263.0) #float(293.15)
    windSpeed10m = float(0.0)
    windDirection10m = float(0.0)
    # assuming this is zenith? probably need Az too?
    main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m)
 
