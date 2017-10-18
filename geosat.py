# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 12:47:03 2017
Geometric Surface Calculation
@author: Danilo
"""

import numpy as np
from quateul import eulerVSquat as quat
from GeometricShader import *


# Geometry Surface Description for 
class SurfaceGeometry:
 
    
    
    # init procedur
    def __init__(self):
        self._q = quat()
    # Global Privat Name
        self._name       = 'surface'
        self._ident      = 1
  
    # Vector / DCM for Orientation
        self._vori       = np.array([0,0,1])
        self._dcmori     = np.array([[1,0,0],
                            [0,1,0],
                            [0,0,1]])
    # Surface Description
        self._alpha      = 1.
        self._epsilon    = 1.
        self._alphaSol   = 1.
        self._epsilonSol = 1.
        self._solararea  = 0.
        self._area       = 2.
        self._check      = 0
        self._efficiency = 0.25
        self._spacer = '_'
        self._shader =     dummyShader()
        self._fullname   = self._name + self._spacer + str(self._ident)

        
    # init simple face with vector, name and id
    def initSimpleSurface(self, v, name, ident):
        self._name       = name
        self._ident      = ident
        self._fullname   = self._name + self._spacer + str(self._ident)
        # Vector / DCM for Orientation
        self._vori       = v / np.linalg.norm(v)
        self._check      = 0
        
        # init simple face with vector, name and id
    def initShadedSurface(self, dcm, name, ident, shader):
        self._name       = name
        self._ident      = ident
        self._fullname   = self._name + self._spacer + str(self._ident)
        # Vector / DCM for Orientation
        self._dcmori     = self._q.correctDcm(dcm)
        self._vori       = np.array([self._dcmori[0,2],self._dcmori[1,2],self._dcmori[2,2]])
        self._check      = 1
        self._shader     = shader
        
      
        # defines Surface Parameters
    def setSurfaceParam(self, alpha, epsilon, area):
        self._alpha       = alpha
        self._epsilon     = epsilon
        self._area        = area
       
         # defines Surface Parameters with Solar
    def setSurfaceParamSolar(self, alpha, epsilon, area, solararea, alphaSolar, epsilonSolar, efficiency):
        self._alpha       = alpha
        self._epsilon     = epsilon
        self._area        = area
        self._solararea   = solararea   
        self._alphaSol    = alphaSolar
        self._epsilonSol  = epsilonSolar
        self._efficiency  = efficiency 
                
        # defines Optical Parameters        
    def setOpticalParam(self, alpha, epsilon):
        self._alpha       = alpha
        self._epsilon     = epsilon
        
        # defines Size
    def setSize(self, area, solararea):
        self._area        = area
        self._solararea   = solararea
     
        # define Solar Optics
    def setOpticalParamSolar(self, alphaSolar, epsilonSolar, efficiency):
        self._alphaSol    = alphaSolar
        self._epsilonSol  = epsilonSolar
        self._efficiency  = efficiency 
        
    def setSpacer(self, spacer):
        self._spacer = spacer
        
        

    
 
    