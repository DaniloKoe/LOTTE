# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:10:22 2017
One Node Solution for Thermal Preprocessing
@author: Danilo Költzsch
"""
import numpy as np
from geosat import SurfaceGeometry as sfg


 # Simple Solver for mid Temperature for Detection of extrem Scenarios
class ThermalOneNodeSolver:
    
       # inits solver by definig geometry
    def __init__(self):
        self._bolzmann         = 5.67 * 10**-8
        self._TKelvin          = 273.15
        self._aluSpeciCapacity = 951
        self._mass             = 1.0
        self._thermalCapacity  = self._mass * self._aluSpeciCapacity
        self._tincr            = 10
        self._count            = 2
        self._sf               = [sfg()]
        self._faceCount        = 1
        
        
       #  redefining surfacegeometry and mass
    def defParam(self, surfaceGeometry, mass):
        self._sf               = surfaceGeometry
        self._mass             = mass
        self._thermalCapacity  = self._mass * self._aluSpeciCapacity
        self._faceCount        = len(surfaceGeometry) 
        
        # calculates middletemperature by init Temperature and previous calculated heatflux and dissipation for all timesteps
    def calcTemp(self, heatflux, dissipation, initTemp):
        dissCount = dissipation.shape
            # objects must have same amount of timesteps
        if heatflux._count == dissCount[0]:
            self._count = heatflux._count
            self._tincr = heatflux._tincr
            self._container   = np.zeros([self._count,6])
            # Temperature
            self._Temp      = np.zeros([self._count])
            self._Temp[0]   = initTemp 
            # Work / Netto / Dissipation / Incomming Radiation / Outgoing Radiation
            self._Wnetto    = np.zeros([self._count])
            self._Wdissi    = np.zeros([self._count])
            self._Win       = np.zeros([self._count])
            self._Wout      = np.zeros([self._count])
            self._WinArea   = np.zeros([self._count, self._faceCount])
            self._WoutArea  = np.zeros([self._count, self._faceCount])
            for i in range(0, (self._count)):
                tempWin = 0
                tempWout = 0
                # calc for every surface of incomming and outgoing radiation
                for j in range(0, self._faceCount):
                    tempWin  = self._Win[i]
                    tempWout = self._Wout[i]
                    self._WinArea[i,j]   = heatflux._heatsum[i,j] * self._sf[j]._alpha * self._sf[j]._area * 10**-6
                    self._WoutArea[i,j]  = self._bolzmann * self._sf[j]._epsilon * self._sf[j]._area * (self._Temp[i] + self._TKelvin)**4 * 10**-6
                    self._Win[i]         = self._WinArea[i,j] + tempWin
                    self._Wout[i]        = self._WoutArea[i,j] + tempWout
                self._Wdissi[i] = dissipation[i]
                self._Wnetto[i] = self._Wdissi[i] + self._Win[i] - self._Wout[i]
                # for the last step
                if i < (self._count - 1):
                    self._Temp[i + 1] = self._Wnetto[i] * self._tincr / self._thermalCapacity + self._Temp[i]
            self._container[:,0] = heatflux._tins[:]
            self._container[:,1] = self._Temp[:]
            self._container[:,2] = self._Wnetto[:]
            self._container[:,3] = self._Wdissi[:]
            self._container[:,4] = self._Win[:]
            self._container[:,5] = self._Wout[:]   
        else:
            print('Heatflux and Dissipation have diffrent amount of timesteps')
            self._container   = np.zeros([self._count,6])
        
                # calculates middletemperature by init Temperature and previous calculated heatflux and dissipation for all timesteps
    def calcTempRed(self, heatflux, dissipation, initTemp):
        dissCount = dissipation.shape
            # objects must have same amount of timesteps
        if heatflux._count == dissCount[0]:
            self._count = heatflux._count
            self._tincr = heatflux._tincr
            self._container   = np.zeros([self._count,6])
            # Temperature
            self._Temp      = np.zeros([self._count])
            self._Temp[0]   = initTemp 
            # Work / Netto / Dissipation / Incomming Radiation / Outgoing Radiation
            self._Wnetto    = np.zeros([self._count])
            self._Wdissi    = np.zeros([self._count])
            self._Win       = np.zeros([self._count])
            self._Wout      = np.zeros([self._count])
            self._WinArea   = np.zeros([self._count, self._faceCount])
            self._WoutArea  = np.zeros([self._count, self._faceCount])
            for i in range(0, (self._count)):
                tempWin = 0
                tempWout = 0
                # calc for every surface of incomming and outgoing radiation
                for j in range(0, self._faceCount):
                    tempWin  = self._Win[i]
                    tempWout = self._Wout[i]
                    self._WinArea[i,j]   = heatflux._redheatsum[i,j] * self._sf[j]._alpha * self._sf[j]._area * 10**-6
                    self._WoutArea[i,j]  = self._bolzmann * self._sf[j]._epsilon * self._sf[j]._area * (self._Temp[i] + self._TKelvin)**4 * 10**-6
                    self._Win[i]         = self._WinArea[i,j] + tempWin
                    self._Wout[i]        = self._WoutArea[i,j] + tempWout
                self._Wdissi[i] = dissipation[i]
                self._Wnetto[i] = self._Wdissi[i] + self._Win[i] - self._Wout[i]
                # for the last step
                if i < (self._count - 1):
                    self._Temp[i + 1] = self._Wnetto[i] * self._tincr / self._thermalCapacity + self._Temp[i]
            self._container[:,0] = heatflux._tins[:]
            self._container[:,1] = self._Temp[:]
            self._container[:,2] = self._Wnetto[:]
            self._container[:,3] = self._Wdissi[:]
            self._container[:,4] = self._Win[:]
            self._container[:,5] = self._Wout[:]   
        else:
            print('Heatflux and Dissipation have diffrent amount of timesteps')
            self._container   = np.zeros([self._count,6])
        
  