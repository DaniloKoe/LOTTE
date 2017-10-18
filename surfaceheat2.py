# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 15:43:46 2017
Heat Flux on Surface Calculation based on STEM
@author: Danilo
"""

from orbit import orbit
#from moi import MoInertia as moi
from sat import satorient as sato
from quateul import eulerVSquat as quat
import numpy as np
from geosat import SurfaceGeometry as sfg

class surfaceheatflux:
   
    
    def __init__(self):
        self._q   = quat()
        self._orb = orbit()
        self._sat = sato()
        self._tins = np.zeros(1)
        self._sf  = [sfg(),sfg(),sfg(),sfg(),sfg(),sfg()]
    # heatflux of the incomming sun (solar spectrum)
        self._heatsun = np.zeros([1,6])
    # heatflux of the incomming reflected sun by earth (albedo) (solar spectrum)
        self._heatalb = np.zeros([1,6])
    # heatflux of the incomming earth infrared (ir spectrum)
        self._heateir = np.zeros([1,6])
    # combined heatfluxes
        self._heatsum = np.zeros([1,6])
    # reduced heatflux of the incomming sun (solar spectrum) (solar to electricity)
        self._redheatsun = np.zeros([1,6])
    # reduced heatflux of the incomming reflected sun by earth (albedo) (solar spectrum) (solar to electricity)
        self._redheatalb = np.zeros([1,6])
    # reduced heatflux of the incomming sun (solar + ir spectrum) (solar to electricity)
        self._redheatsum = np.zeros([1,6])
        self._anglesun = np.zeros([1,6])
        self._angleearth = np.zeros([1,6])
    
    # power absorbed from different faces 
        self._powersun = np.zeros([1, 6])
        self._poweralb = np.zeros([1, 6])
        self._powereir = np.zeros([1, 6])
        self._powersum = np.zeros([1, 6])
        self._sumpowersol  =  np.zeros([1])
        self._sumpoweralb  =  np.zeros([1])
        self._sumpowereir  =  np.zeros([1])
        self._sumpowersum  =  np.zeros([1])
        self._sumpowersolarrad  =  np.zeros([1])
    # electrical power generated from solar generators
        self._solarenergy  =  np.zeros([1, 6])
    # consumed solar energy (used from solar generators)
        self._powerelec =  np.zeros([1, 6])
    # reduced sum power
        self._redpowersum = np.zeros([1, 6])
        self._sumsolarenergy  =  np.zeros([1])
    
    # sum power solar for electricity and for heat
        self._sumpowerelec  =  np.zeros([1])
        self._sumpowerheat  =  np.zeros([1])
    # ecosump for solar cell reduction
        self._econsump  =  np.zeros([1])
    # Marker if reduced heatflux is performed
        self._redmark  = 0

    
  #      self._vector   = np.zeros([1,6])
        self._nface = 6
        self._count = 2
        self._tincr = 10
    # Possible Cases (hoz, medium, cold)
        self._caseSun = 'medium'
        self._caseAlb = 'medium'
        self._caseEIr = 'medium'
    # Hot Cases in W/m^2 / W/m^2 and as Faktor
        self._SolarHC  = 1414.
        self._EIrHC    = 232.
        self._AlbedoHC = 0.23
    # Medium Cases in W/m^2 / W/m^2 and as Faktor
        self._SolarMC  = 1367.
        self._EIrMC    = 228.
        self._AlbedoMC = 0.215
    # Cold Cases in W/m^2 / W/m^2 and as Faktor
        self._SolarCC  = 1322.
        self._EIrCC    = 224.
        self._AlbedoCC = 0.2
    # Choice Standard Medium Case
        self._SolarC   =     self._SolarMC 
        self._EIrC     =     self._EIrMC
        self._AlbedoC  =     self._AlbedoMC
        #print('Heatflux Object Inited')
    
    # Definition of Heatflux Environment by Numbers in W/m^2 / W/m^2 and as Faktor
    def setEnvironmentByNumber(self, SolarCt, EIrCt,  AlbedoCt):
        self._SolarC   = SolarCt
        self._EIrC     = EIrCt
        self._AlbedoC  = AlbedoCt
    
    # Definition of Heatflux Environment by Possible Cases (hoz, medium, cold)
    def setEnvironmentByCase(self, SolarCase, EIrCase,  AlbedoCase):
        # define Solar Case
        if SolarCase == 'hot':
            self._SolarC   = self._SolarHC
        elif SolarCase == 'medium':
            self._SolarC   = self._SolarMC
        elif SolarCase == 'cold':
            self._SolarC   = self._SolarCC
        else:
            self._SolarC   = self._SolarMC
            print('Wrong Solar Case defined, medium Case chosen')
        # define Earth Infrared Case    
        if EIrCase == 'hot':
            self._EIrC   = self._EIrHC
        elif EIrCase == 'medium':
            self._EIrC   = self._EIrMC
        elif EIrCase == 'cold':
            self._EIrC   = self._EIrCC
        else:
            self._EIrC   = self._EIrMC
            print('Wrong Eart Infrared Case defined, medium Case chosen')
        # define Albedo Case    
        if AlbedoCase == 'hot':
            self._AlbedoC   = self._AlbedoHC
        elif AlbedoCase == 'medium':
            self._AlbedoC   = self._AlbedoMC
        elif AlbedoCase == 'cold':
            self._AlbedoC   = self._AlbedoCC
        else:
            self._AlbedoC   = self._AlbedoMC
            print('Wrong Albedo Case defined, medium Case chosen')    

    # Surface Heatflux Calc
    def surfaceHeatflux(self, orbit, orient, geosat):
        # Surface Heatflux Implements by defined Data
        self._orb = orbit
        earthr    = self._orb._rearth
        self._sat = orient
        self._sf  = geosat
        self._count = self._sat._count
        self._tincr = self._orb._tincr
        self._nface = len(self._sf)
        # Set local V
        tempsatv    = np.zeros(3)
        tempsunsatv = np.zeros(3)
        #tempsunv    = np.zeros(3)
        tempvface   = np.zeros(3)
        temporientv = np.zeros(10)
        tempvface   = np.zeros(3)
        tempdcmface = np.zeros([3,3])
        self._tins  = np.zeros(self._count)
        # heatfluxes
        self._heatsun = np.zeros([self._count, self._nface])
        self._heatalb = np.zeros([self._count, self._nface])
        self._heateir = np.zeros([self._count, self._nface])
        # combined heatfluxes
        self._heatsum = np.zeros([self._count, self._nface])
        # angle between normal of area and incomming radiation
        self._angleearth =  np.zeros([self._count, self._nface])
        self._anglesun   =  np.zeros([self._count, self._nface])
        
        # reduced heat flux
        self._redheatsun = np.zeros([self._count, self._nface])
        self._redheatalb = np.zeros([self._count, self._nface])
        self._redheatsum = np.zeros([self._count, self._nface])
        # power absorbed from areas
        self._powersun = np.zeros([self._count, self._nface])
        self._poweralb = np.zeros([self._count, self._nface])
        self._powereir = np.zeros([self._count, self._nface])
        self._powersum = np.zeros([self._count, self._nface])
        self._powerelec =  np.zeros([self._count, self._nface])
        self._redpowersum = np.zeros([self._count, self._nface])
        # power absorbed from all areas
        self._sumpowersol  =  np.zeros([self._count])
        self._sumpoweralb  =  np.zeros([self._count])
        self._sumpowereir  =  np.zeros([self._count])
        self._sumpowersum  =  np.zeros([self._count])
        
        # electrical power generated from solar generators
        self._solarenergy  =  np.zeros([self._count, self._nface])
        self._sumsolarenergy  =  np.zeros([self._count])
        
        # sum power solar for electricity and for heat
        self._sumpowerelec  =  np.zeros([self._count])
        self._sumpowerheat  =  np.zeros([self._count])
        
        # time step
        for i in range(0, self._count):
            # reads local fields for timpe step
            self._tins[i]  = self._orb._tins[i]
            tempsatv = self._orb._satv[i,:]
            tempsatr = np.linalg.norm(tempsatv)
            tempsunsatv = self._orb._sunsatv[i,:]            
            tempillu = self._orb._light[i]
            #tempsunv = self._orb._sunv[i,:]            
            temporientv = self._sat._orientV[i,:]
            temporidcm  = np.array([[temporientv[1], temporientv[4], temporientv[7]],
                                    [temporientv[2], temporientv[5], temporientv[8]],
                                    [temporientv[3], temporientv[6], temporientv[9]]])
            tempshadearth = 1
            tempshadesun  = 1
            # surfaces
            for j in range(0, self._nface):
                # reads and rotats surface in sat orientation
                sfgtemp = self._sf[j]
                # for unshaded faces
                
                # Local Orbit Vector for each time step
       # self._vector
                
                if sfgtemp._check == 0:
                    tempvface[:] = sfgtemp._vori[:]
                    tempvface2   = self._q.DcmRotOfV(temporidcm, tempvface)
                    # calc for solar radiantion on surface, negative direction for angle calc
                    cosanglesun   =  self._q.getRadiCosByV(-tempsunsatv, tempvface2)
                    beta =  self._q.getAngleRadByV(-tempsatv, tempvface2)
                    self._angleearth[i,j]   = beta  
                    h = tempsatr/ earthr
                    # is sphere by phase cutted
                    casesphere = np.pi / 2 - np.math.asin(1/h)
                    # criteria for end of calc of viewfactor
                    if h < 1.72:
                        endview = 12.018*h**2-28.35*h+18.96
                    else:
                        endview = 2.41
                        
                    if beta < casesphere:
                        viewearth = np.math.cos(beta)/h**2
                    elif beta < endview:
                        x = np.sqrt(h**2 - 1)
                        y = - x / np.math.tan(beta)
                        z = np.math.sin(beta * np.sqrt(1 - y**2))
                        viewearth = 1 / ( np.pi * h**2) * (np.math.cos(beta) * np.math.acos(y) - x * z) + 1 / np.pi * np.math.atan(z / x)
                    else:
                        viewearth = 0.

                    self._heatsun[i,j] = cosanglesun * tempillu * self._SolarC
                    self._heateir[i,j] = viewearth *  self._EIrC 
                    self._heatalb[i,j] = viewearth * tempillu * self._SolarC * self._AlbedoC
                    self._heatsum[i,j] = self._heatsun[i,j] + self._heateir[i,j] + self._heatalb[i,j]
                    if cosanglesun < 1:                        
                        self._anglesun[i,j]  = np.arccos(cosanglesun) * 180 / np.pi
                    else:
                        self._anglesun[i,j]  = 0.0    
                    # heat calc defined by radiation and geometric class
                    self._powersun[i,j] =  self._heatsun[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._powereir[i,j] =  self._heateir[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._poweralb[i,j] =  self._heatalb[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._powersum[i,j] =  self._powersun[i,j] + self._powereir[i,j] + self._poweralb[i,j]
                    self._sumpowersol[i] = self._sumpowersol[i] + self._powersun[i,j] 
                    self._sumpowereir[i] = self._sumpowereir[i] + self._powereir[i,j] 
                    self._sumpoweralb[i] = self._sumpoweralb[i] + self._poweralb[i,j] 
                    self._sumpowersum[i] = self._sumpowersum[i] + self._powersum[i,j]
                    # calc of electric energy absorbed
                    tempsolarheatflux   = self._heatsun[i, j] + self._heatalb[i, j] 
                    self._solarenergy[i,j] = tempsolarheatflux * sfgtemp._efficiency * sfgtemp._solararea / 10.**6
                    self._sumsolarenergy[i] = self._sumsolarenergy[i]  + self._solarenergy[i,j] 
                    
                elif sfgtemp._check == 1:
                    tempdcmface[:,:] = sfgtemp._dcmori[:,:]
                    tempdcmface2 = self._q.DcmRotOfDcm(temporidcm, tempdcmface)
                    # calc for solar radiantion on surface, negative direction for angle calc
                    cosanglesun   =  self._q.getRadiCosByVZofDcm(-tempsunsatv, tempdcmface2)
                    beta =  self._q.getAngleRadByV(-tempsatv, tempvface2)
                    self._angleearth[i,j]   = beta  
                    h = tempsatr/ earthr
                    # is sphere by phase cutted
                    casesphere = np.pi / 2 - np.math.asin(1/h)
                    # criteria for end of calc of viewfactor
                    
                    if h < 1.72:
                        endview = 12.018*h**2-28.35*h+18.96
                    else:
                        endview = 2.41
                        
                    if beta < casesphere:
                        viewearth = np.math.cos(beta)/h**2
                    elif beta < endview:
                        x = np.sqrt(h**2 - 1)
                        y = - x / np.math.tan(beta)
                        z = np.math.sin(beta * np.sqrt(1 - y**2))
                        viewearth = 1 / ( np.pi * h**2) * (np.math.cos(beta) * np.math.acos(y) - x * z) + 1 / np.pi * np.math.atan(z / x)
                    else:
                        viewearth = 0.
                        
                    angle1sun, angle2sun  =  self._q.getShaderAngleByVZofDcm(-tempsunsatv, tempdcmface2)
                    angle1earth, angle2earth  =  self._q.getShaderAngleByVZofDcm(-tempsatv, tempdcmface2)
                    factsun   = sfgtemp._shader.calcShader(angle1sun, angle2sun)
                    factearth = sfgtemp._shader.calcShader(angle1earth, angle2earth)
                    # here has to be shader
                    # radiation calc, defined by Jakob Kühn
                    self._heatsun[i,j] = factsun * cosanglesun * tempillu * self._SolarC
                    self._heateir[i,j] = factearth * viewearth *  self._EIrC
                    self._heatalb[i,j] = factearth * viewearth * tempillu * self._SolarC * self._AlbedoC
                    self._heatsum[i,j] = self._heatsun[i,j] + self._heateir[i,j] + self._heatalb[i,j]
                    if cosanglesun < 1:                        
                        self._anglesun[i,j]  = np.arccos(cosanglesun) * 180 / np.pi
                    else:
                        self._anglesun[i,j]  = 0.0 
                    # heat calc defined by radiation and geometric class
                    self._powersun[i,j] =  self._heatsun[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._powereir[i,j] =  self._heateir[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._poweralb[i,j] =  self._heatalb[i,j] *  sfgtemp._alpha * sfgtemp._area / 10.**6
                    self._powersum[i,j] =  self._powersun[i,j] + self._powereir[i,j] + self._poweralb[i,j]
                    self._sumpowersol[i] = self._sumpowersol[i] + self._powersun[i,j] 
                    self._sumpowereir[i] = self._sumpowereir[i] + self._powereir[i,j] 
                    self._sumpoweralb[i] = self._sumpoweralb[i] + self._poweralb[i,j] 
                    self._sumpowersum[i] = self._sumpowersum[i] + self._powersum[i,j]
                    # calc of electric energy absorbed
                    tempsolarheatflux   = self._heatsun[i, j] + self._heatalb[i, j] 
                    self._solarenergy[i,j] = tempsolarheatflux * sfgtemp._efficiency * sfgtemp._solararea / 10.**6
                    self._sumsolarenergy[i] = self._sumsolarenergy[i]  + self._solarenergy[i,j] 
                    
                else:
                    print('Error Wrong Shader Case Defined')
        print('Heat Fluxes Defined for all Surfaces and Steps')

        
        
        # calculates solar energy for all faces
    def getSolarenergy(self, timestep):
        solEnergy = np.zeros(self._nface)
        # reads all previous solar energies for one time step
        solEnergy[:] = self._solarenergy[timestep,:]
        return solEnergy
    
            # calculates solar energy for all faces for all steps
    def getAllSolarenergy(self):
        # reads all previous solar energies for one time step
        return self._solarenergy

    
    
         # calculates solar energy for all faces
    def sumSolarenergy(self, timestep):
        solEnergy = np.zeros(self._nface)
        # reads all previous solar energies for one time step
        solEnergy[:] = self._solarenergy[timestep,:]
        return solEnergy
            
        
     
        # Surface Heatflux reducment by transformed solar radiation
    def reduceHeatflux(self, solarcellload):
        self._econsump = np.zeros([self._count])
        loccount = np.shape(solarcellload)
        # checks if field lengths matches
        if self._count ==  loccount[0]:
            self._redmark = 1
            self._econsump[:] = solarcellload[:]
            for i in range(0, self._count):
                penaltyfaktor = 1
                tempsolenergy = self._sumsolarenergy[i]
                # calculation of reduction factor relating energy tranformation into electricity to total absorbed solar energy
                # if no solar energy is incomming
                if tempsolenergy > 0:
                    penaltyloadfaktor = (solarcellload[i]) / tempsolenergy 
                else:
                    penaltyloadfaktor = 1
                # reduction of heatfluxes for every surface
                for j in range(0, self._nface):
                    tempsumsol            = self._powersun[i,j] + self._poweralb[i,j]
                    # relates generated solar energy with the incomming solar radiation
                    self._powerelec[i,j]  = penaltyloadfaktor * self._solarenergy[i,j]
                    if tempsumsol > 0 and tempsumsol > self._powerelec[i,j]:
                        penaltyfaktor     = (tempsumsol - self._powerelec[i,j]) / tempsumsol
                    else:
                        penaltyfaktor     = 1
                    self._redheatsun[i,j] = penaltyfaktor * self._heatsun[i,j] 
                    self._redheatalb[i,j] = penaltyfaktor * self._heatalb[i,j] 
                    self._redheatsum[i,j] = self._redheatsun[i,j] + self._redheatalb[i,j] + self._heateir[i,j]
                    self._redpowersum[i,j] = self._powersum[i,j] - self._powerelec[i,j]
                # sum power solar for electricity and for heat
                self._sumpowerelec[i]  =  solarcellload[i]
                self._sumpowerheat[i]  =  self._sumpowersum[i] - solarcellload[i]
            print('Solar Heatfluxes reduced respecting solar energy production')
        else:
            print('ERROR: Shape does not match')
                
                
                
            # Surface Heatflux reducment by transformed solar radiation [case: 0 - standard heatflux [W/m^2]; 1 - reduced heatflux  [W/m^2]; 2 - standard power [W]; 3 - reduced power [W]; 4 - standard albedo respecting heatflux [W/m^2]; 5 - reduced albedo respecting heatflux [W/m^2]]
    def getHeatflux(self, case):
        # case:  - reduced heatflux  [W/m^2]
        if case == 1 and self._redmark == 1:
            return self._redheatsum
        # case: 0 - standard heatflux [W/m^2]
        elif case == 1 and self._redmark == 0:
            print('ERROR: No reduceHeatflux performed')
            return self._redheatsum
        # case: 2 - standard power [W]
        elif case == 2:
            return self._powersum
        # case: 3 - reduced power [W]
        elif case == 3 and self._redmark == 1:
            return self._redpowersum
        # case: 2 - standard power [W]
        elif case == 3 and self._redmark == 0:
            print('ERROR: No reduceHeatflux performed')
            return self._powersum
        # case: 4 - standard absorbing factor (alpha) respecting heatflux [W/m^2]
        elif case == 4:
            albheatflux = np.zeros([self._count, self._nface])
            for i in range (0, self._nface):
                albheatflux[:,i] = self._sf[i]._alpha * self._heatsum[:,i]
            return albheatflux
        # case:  5 - reduced absorbing factor (alpha) respecting heatflux [W/m^2]
        elif case == 5 and self._redmark == 1:
            albheatflux = np.zeros([self._count, self._nface])
            for i in range (0, self._nface):
                albheatflux[:,i] = self._sf[i]._alpha * self._redheatsum[:,i]
            return albheatflux
        # case: 4 - standard absorbing factor (alpha) respecting heatflux [W/m^2]
        elif case == 5 and self._redmark == 0:
            albheatflux = np.zeros([self._count, self._nface])
            for i in range (0, self._nface):
                albheatflux[:,i] = self._sf[i]._alpha * self._heatsum[:,i]
            print('ERROR: No reduceHeatflux performed')
            return albheatflux    
        # case: 0 - standard heatflux [W/m^2]
        else:
            return self._heatsum
                
        
                    
                    
                    
                    
                    
    
        
        
        
        
    
