# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:47:19 2017
Moment of Inertia Caclulation
@author: Danilo Költzsch
"""

import numpy as np
from scipy.linalg import eig

class MoInertia:

     
       
    
    
       def __init__(self):
                   # matrix values in g * square metre therfore / 1000 to be in kg/m^2
           self._jtensor = np.matrix([[1.8,   0,   0],
									  [0,   1.8,   0],
									  [0,     0, 2.2]]) / 1000
           self._evm = np.matrix([[1,       0,       0],
                        [0,       1,       0],
                        [0,       0,       1]]) 
         # matrix values in g * square metre therfore / 1000 to be in kg/m^2
           self._princ = np.matrix([[[1.8,   0,    0],
									  [0,   1.8,   0],
									  [0,     0, 2.2]]) / 1000
           self._calc = 'false'
           #print('MoI Object initialised')
       
       # calcs the principal axis of moi and their values 
       def calceig(self):
           ew,ev = eig(self._jtensor)
           self._evm = np.matrix([[np.sign(ev[0,0]) * ev[0,0], np.sign(ev[1,1]) * ev[0,1], np.sign(ev[2,2]) * ev[0,2]],
                                 [ np.sign(ev[0,0]) * ev[1,0], np.sign(ev[1,1]) * ev[1,1], np.sign(ev[2,2]) * ev[1,2]],
                                 [ np.sign(ev[0,0]) * ev[2,0], np.sign(ev[1,1]) * ev[2,1], np.sign(ev[2,2]) * ev[2,2]]])
           self._ewm = np.matrix([[ew[0],      0,        0],
                                  [       0, ew[1],      0],
                                  [       0,      0, ew[2]]])
           self._calc = 'true' 
           print("Moment of Inertia calculated")
           return self._ewm, self._evm
       
       # dir calcs the principal axis of moi and their values 
       def dirCalcEig(self, jtensor):
           princ,ev = eig(jtensor)           
           self._evm = np.matrix([[np.sign(ev[0,0]) * ev[0,0], np.sign(ev[1,1]) * ev[0,1], np.sign(ev[2,2]) * ev[0,2]],
                                 [ np.sign(ev[0,0]) * ev[1,0], np.sign(ev[1,1]) * ev[1,1], np.sign(ev[2,2]) * ev[1,2]],
                                 [ np.sign(ev[0,0]) * ev[2,0], np.sign(ev[1,1]) * ev[2,1], np.sign(ev[2,2]) * ev[2,2]]])
           self._princ = np.matrix([[princ[0],      0,        0],
                                  [       0, princ[1],      0],
                                  [       0,      0, princ[2]]])
           self._calc = 'true' 
           print("Moment of Inertia calculated")
           return self._princ , self._evm 
        
          
       # defines new tensor, needs 3x3 matrix
       def newJtensor(self, newjtens):
           self._jtensor = np.matrix([[newjtens[0,0],newjtens[0,1],newjtens[0,2]],
                                      [newjtens[1,0],newjtens[1,1],newjtens[1,2]],
                                      [newjtens[2,0],newjtens[2,1],newjtens[2,2]]])
    
       # defines new tensor, needs components
       def newJcomp(self, jxx, jyy, jzz, jxy, jxz, jyz):
           self._jtensor = np.matrix([[jxx, -jxy, -jxz],
                                     [ -jxy, jyy, -jyz],
                                     [ -jxz, -jyz, jzz]])
    
       #returns principal Moment of Inertia  relativ to coordinate system and center of gravity
       def getMoI(self):
           self.calceig()
           return self._princ, self._evm
               
       #returns MoI Tensor relativ to coordinate system and center of gravity
       def getJ(self):
           return self._jtensor       