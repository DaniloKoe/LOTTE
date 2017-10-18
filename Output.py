# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 14:21:55 2017
Generates Output for further exterior Caculation
@author: Danilo Költzsch
"""

#from orbit import orbit
#from moi import MoInertia as moi
#from sat import satorient as sato
#from quateul import eulerVSquat as quat

from geosat import SurfaceGeometry as sfg
import numpy as np
from surfaceheat import surfaceheatflux as sfh
from openLibrary import *
import os
import pandas as pd
import time



class OutputGenerator:
   
    
        # inits Outp8ut generator with objects of heatflux and esystem library object
    def __init__(self, heatflux, esystem, note):
        self._nameheatfluxf = 'sf_'
        self._powerheaderadd = ' in W'
        self._heatfheaderadd = ' in W/m^2'
        self._dissheaderadd  = ' in W'
        self._heatheaderadd  = ' in W/m^2'
        self._reducheaderadd = ' reduced'
        self._albedoadd      = ' albedo inlc.'
        self.ANSYSext        = '.txt'
        self.ANSYSl          = 'loadstep'
        self.ANSYS2          = 'timestep'
        self.ANSYS3          = 'power_'
        self.ANSYS4          = 'heatf_'
        self.ANSdissadd      = self.ANSYS3
        self.ANSheatadd      = ''
        self.ANSdissicase    = 0
        self.ANSheatfcase    = 1
        self._esystem  = esystem
        self._heatflux = heatflux
        self._count = self._heatflux._count
        self._nface = self._heatflux._nface
        self._disname = self._esystem.getName()
        self._ndis    = len(self._disname)
        self._diss    = np.zeros([self._count, self._ndis])
        self._sfheat  = np.zeros([self._count, self._nface])
        self._ANSdiss    = np.zeros([self._count, self._ndis])
        self._ANSsfheat  = np.zeros([self._count, self._nface])
        self._tincr   = self._heatflux._tincr
        self._tins    = np.zeros([self._count])
        self._tins[:] = self._heatflux._tins[:]
        # Standard
        self._sfheat[:,:]  = self._heatflux._heatsum[:,:]
        self._diss         = self._esystem.getDissipation()
        self._note         =  note + ' ' + time.strftime("%H-%M %d.%m.%Y")
           


                       
        # writes heatflux with case [ [0] standard heatflux [W/m^2], [1] reduced heatflux [W/m^2],  [2] standard power (related to defined surface [W],  [3] reduced power (related to defined surface [W]; [4] standard heatflux albedo respecting [W/m^2], [5] reduced heatflux albedo respecting [W/m^2]]
    def heatfluxCase(self, case):
         #case:   [0] standard heatflux [W/m^2]
        self._sfheat[:,:]  = self._heatflux.getHeatflux(case)
        self.ANSheatadd = ''
        self.ANSheatfcase   = 1
         #case:   [1] reduced heatflux [W/m^2]
        if case == 1:
            self._heatheaderadd  = self._reducheaderadd + self._heatfheaderadd
            self.ANSheatfcase   = 1
         #case:   [2] standard power (related to defined surface) [W]
        elif case == 2:
            self._heatheaderadd  = '' + self._powerheaderadd
            self.ANSheatadd      = self.ANSYS3
            self.ANSheatfcase   = 0
         #case:   [3] reduced power (related to defined surface) [W]
        elif case == 3:
            self._heatheaderadd  = self._reducheaderadd + self._powerheaderadd
            self.ANSheatadd      = self.ANSYS3
            self.ANSheatfcase    = 0
         #case:   [4] standard heatflux absorbing factor (alpha) respecting  (related to defined surface)[W/m^2]
        elif case == 4:
            self._heatheaderadd  = self._albedoadd + self._heatfheaderadd
            self.ANSheatfcase    = 1
         #case:   [5] reduced heatflux absorbing factor (alpha) respecting  (related to defined surface)[W/m^2]
        elif case == 5:
            self._heatheaderadd  = self._albedoadd + self._reducheaderadd + self._heatfheaderadd
            self.ANSheatfcase    = 1
        

            
        # returnes Dissipation defined as Power in [W] [0] or as heatflux related to the surface size [W/m^2]
    def dissipationCase(self, case):
        if case == 1.:
            #case:  [1] dissipation defined as heatflux  (related to defined surface)[W/m^2]
            self._diss  = self._esystem.getHeatflux()
            self._dissheaderadd  = '' + self._heatfheaderadd
            self.ANSdissadd     = self.ANSYS4 
            self.ANSdisscase    = 1
        else:
            #case:   [0] standard dissipation  [W]
            self._diss  = self._esystem.getDissipation()
            self._dissheaderadd  = '' + self._powerheaderadd
            self.ANSdissadd     = self.ANSYS3
            self.ANSdisscase    = 0
            
            
        # reduction for all relevant components for less timesteps    
    def reduction(self, stepnew):
        self._diss   = ReduceOutput.reduceMultiArray(self._diss, stepnew)
        self._sfheat = ReduceOutput.reduceMultiArray(self._sfheat, stepnew)
        self._tins   = ReduceOutput.reduceTimestep(self._tins, stepnew)
        self._tincr  = self._tins[0] - self._tins[1]
        self._count  = stepnew + 1
        
        
        
    def writeXLSX(self, name, file_path):
        cols = ['time [s]']
        columns = (self._nface + self._ndis + 1)
        plotmatrix = np.zeros([self._count, columns])
        plotmatrix[:,0] = self._tins[:]
        #set file for save
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        for i in range(0, self._nface):
            headplot = self._heatflux._sf[i]._fullname + self._heatheaderadd
            cols.append(headplot)
            j = i + 1
            plotmatrix[:,j] = self._sfheat[:,i]
        for i in range(0, self._ndis):
            headplot =  self._disname[i] + self._dissheaderadd
            cols.append(headplot)
            j = i + 1 + self._nface
            plotmatrix[:,j] = self._diss[:,i]
        df = pd.DataFrame(data=plotmatrix, columns=cols)
        writer = pd.ExcelWriter(file_path+name+'.xlsx', engine='xlsxwriter')
        df.to_excel(writer, sheet_name=self._note)
        writer.save()
        print('XLSX File is exported in ' +  name + '.xlsx')
    
    
    def writeCSV(self, name, file_path, seperator):
        cols = ['time [s]']
        columns = (self._nface + self._ndis + 1)
        plotmatrix = np.zeros([self._count, columns])
        plotmatrix[:,0] = self._tins[:]
        #set file for save
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        for i in range(0, self._nface):
            headplot = self._heatflux._sf[i]._fullname + self._heatheaderadd
            cols.append(headplot)
            j = i + 1
            plotmatrix[:,j] = self._sfheat[:,i]
        for i in range(0, self._ndis):
            headplot =  self._disname[i] + self._dissheaderadd
            cols.append(headplot)
            j = i + 1 + self._nface
            plotmatrix[:,j] = self._diss[:,i]
        df = pd.DataFrame(data=plotmatrix, columns=cols)        
        # encoding='ansi' in 3.5 not available anymore
        df.to_csv(file_path+name+'.csv', index=True, header=True, encoding='mbcs',  sep=seperator)
        print('CSV File is exported in ' +  name + '.csv')
    
    
    # creates an ANSYS readable file for further calculation
    def writeANSYS(self, name, file_path):
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        if self.ANSdisscase == 1:
            self._ANSdiss = self._diss / 1000
        else:
            self._ANSdiss = self._diss
            
        if self.ANSheatfcase == 1:
            self._ANSsfheat = self._sfheat / 1000
        else:
            self._ANSsfheat = self._sfheat
            
        file = open(file_path+name+self.ANSYSext ,'w') 
        # needed ANSYS Header for interior DIM command (array size)
        file.write('! '+ self._note + '\n')
        file.write(self.ANSYSl + '=' + str(self._count) + '\n')
        file.write('/eof' + '\n')
        file.write('' + '\n')
        file.write('! ' +  name + self.ANSYSext + '\n')
        file.write('' + '\n')
        # number of loadsteps wich shall be solved
        file.write(self.ANSYSl + '=' + str(self._count) + '\n')
        file.write('' + '\n')
        # setting up time array
        file.write('! timesteps in sec:' + '\n')
        for i in range (0, self._count):
            file.write(self.ANSYS2+'('+str(i+1)+')=' + str(self._tins[i]) + '\n') 
        file.write('' + '\n')
        # radiation heat fluxes
        file.write('! heatfluxes on outer surfaces '+ self._heatheaderadd + ' :' + '\n')
        for i in range (0, self._count):
            for j in range(0, self._nface):
                file.write(self.ANSheatadd + self._heatflux._sf[j]._fullname+'('+str(i+1)+')=' + str(self._ANSsfheat[i,j]) + '\n')
        file.write('' + '\n')
        # heat disstribution
        file.write('! inner heat distribution on defined components '+ self._dissheaderadd + ' :' + '\n')
        for i in range (0, self._count):
            for j in range(0, self._ndis):
                file.write(self.ANSdissadd + self._disname[j]+'('+str(i+1)+')=' + str(self._ANSdiss[i,j]) + '\n')
        file.close()
        print('ANSYS Import File is exported in ' +  name + '.txt')
        
        
    

    # reduce for output of arrays
class ReduceOutput:    
    
    # reduces number of fields defined by faktor by building the mean of amount [needs no init]
    def reduceVektor(vektor, stepnew):
        shapev = np.shape(vektor)
        shape  = shapev[0]
        stepnew = int(stepnew)
        step  = stepnew + 1
        mean  = np.zeros(step)
        faktor = stepnew / shape
        if faktor > 0 and faktor <= 1: 
            div      = shape / stepnew
            t = 0
            mean  = np.zeros([step])
            redfac = 1.
            for i in range (0, step):
                if t < shape:
                    j = 1.
                    # dynamic field length for chord - trapezoid mean procedure (Trapez Sehnen Verfahren)
                    if i == 0 or i == stepnew:
                        # reduction factor for half segments (first and last element)
                        fieldle = div / 2.
                    else:
                        fieldle = div
                    while (fieldle - j)  >= 0 and t < shape:
                        mean[i] = mean[i] + redfac * vektor[t]
                        j = j + redfac
                        t = t + 1
                        redfac = 1.
                        # for parted segments
                    if (fieldle - j) > -1  and t < shape:
                        redfac = (j - fieldle)
                        mean[i] = mean[i] + (1 - redfac) * vektor[t]
                    # calculation of segment number
                    mean[i] = mean[i] / fieldle  
        else:
            print('false amount of steps defined')
        return mean
    
    
     # reduces number of fields for a multiple array
    def reduceMultiArray(multiarray, stepnew):
        shapev = np.shape(multiarray)
        shapex  = shapev[1]
        shapey  = shapev[0]
        mean = np.zeros([shapey, shapex])
        step  = stepnew + 1
        faktor = stepnew / shapey
        if faktor > 0 and faktor <= 1: 
            mean  = np.zeros([step, shapex])
            for i in range (0, shapex):
                mean[:,i] = ReduceOutput.reduceVektor(multiarray[:,i], stepnew)
        else:
            print('false amount of steps defined')
        return mean
    
    # reduces the timesteps for timevektor
    def reduceTimestep(timevektor, stepnew):
        shapev = np.shape(timevektor)
        shape  = shapev[0]
        timenew = np.zeros([shape])
        step  = stepnew + 1
        faktor = stepnew / shape
        if faktor > 0 and faktor <= 1: 
            div      = timevektor[-1] / stepnew
            timenew  = np.zeros([step])
            for i in range(0, step):
                timenew[i] = div * i 
            timenew[-1] = timevektor[-1]
        else:
            print('false amount of steps defined')
        return timenew
    
            
        
    
    