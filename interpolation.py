import sys
import numpy as np
from  pycrtm import pycrtm
class profileInterpolate:
    def __init__(self, pOutput, pInput, items):
        """
        Intialize with stuff you want to interpolate 
        pOutput [nlev/nlay] --pressure output vertical grid 
        pInput [nprofiles,nlev/nlay] -- pressure input grid 
        items (dictonary of arrays) [nprofiles,nlev/nlay] -- items on input grid
        """
        self.pOutput = pOutput
        self.pInput = pInput 
        self.items = items
        self.nv = pOutput.shape[0]
        self.nprof = pInput.shape[0]
        self.pOutputGrid = np.zeros([self.nprof,self.nv])
        self.itemsInterp = {}
        for i in list(self.items.keys()):
            if(self.items[i].ndim<3):
                self.itemsInterp[i] = np.zeros([self.nprof,self.nv])
            else:
                self.itemsInterp[i] = np.zeros([self.nprof,self.items[i].shape[1], self.nv])
    def logLinear(self, x, xo, yo):
        """
        Do a log-linear interpolation.
        """
        logX =  np.log(x)
        logXo = np.log(xo)
        logYo = np.log(yo)
        return np.exp(np.interp(logX, logXo, logYo))

    def averageMidpoint(self, po, pi, itemIn):
        """
        Do layer averaging by taking the average between levels.
        """
        itemOut = np.zeros(itemIn.shape[0]-1)
        for n in list( range(1,pi.shape[0]) ):
            itemOut[n-1] = 0.5*(itemIn[n-1]+itemIn[n])
        return itemOut

    def crtmInterpWrap(self, Pout, Pin, Xin):
        """
        Use LayerAvg function from CRTM to Interpolate User pressure Layers to Coefficient levels. 
        """
        Xout = pycrtm.applyavg(np.log(Pout),np.log(Pin),Xin)
        return Xout
     
    def interp(self,pOutput,pInput,itemInput, method='crtm-wrap'):
        """
        Do an interpolation on a profile using selected method.
        """
        if method == 'log-linear':
            return self.logLinear(pOutput,pInput, itemInput)
        elif method == 'average':
            return self.averageMidpoint(pOutput, pInput, itemInput)
        elif method == 'crtm-wrap':
            return  self.crtmInterpWrap(pOutput, pInput, itemInput) 
        else: sys.exit("Unknown Interpolation method")  
    def interpProfiles(self,method='crtm-wrap'):
        """
        Apply interpolation to all provided profiles.
        """
        for i in list( self.items.keys() ):
            for ii in list(range(self.nprof)):
                if(self.items[i].ndim>2): # jacobian with profile, channels, nlevels
                    for k in list(range(self.items[i].shape[1])):
                        self.itemsInterp[i][ii,k,:] = self.interp(self.pOutput, self.pInput[ii,:], self.items[i][ii,k,:], method = method)
                else:
                    self.itemsInterp[i][ii,:] = self.interp(self.pOutput, self.pInput[ii,:], self.items[i][ii,:], method = method)

        # output grid
        if(method=='average'):
            for i in list(range(self.nprof)):
                self.pOutputGrid[i,:] = self.interp(self.pOutput, self.pInput[i,:], self.pInput[i,:], method='average')
        else:
            for i in list(range(self.nprof)):
                self.pOutputGrid[i,:] = self.pOutput     
    def get(self):
        return self.pOutput, self.itemsInterp 
