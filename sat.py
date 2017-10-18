# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 19:36:41 2017
Satellite Orientation Calculation
@author: Danilo
"""

import numpy as np
from orbit import orbit as orb
from moi import MoInertia as moi
from quateul import eulerVSquat as quat

class satorient:
  
    
     # init procedur
    def __init__(self):
       self._jtens = moi()
       self._orbit = orb()
       
       self._q = quat()
       self._count = 2
    # defines MoI Values
       self._J = np.matrix([[1,0,0],
                            [0,2,0],
                            [0,0,3]])
       self._Jprinc = np.matrix([[1,0,0],
                                 [0,2,0],
                                 [0,0,3]])
       self._Jevm = np.matrix([[1,0,0],
                               [0,1,0],
                               [0,0,1]])
    # defines the systems flight direction in local sat system: flight | normal | nadir vector
       self._rotationM =     np.matrix([[1,0,0],
                                        [0,1,0],
                                        [0,0,1]])
       self._startVnadir =   np.array([ self._orbit._satv[0,0],   self._orbit._satv[0,1],   self._orbit._satv[0,2]]) * -1.
       self._startVflight =  np.array([ self._orbit._satvel[0,0], self._orbit._satvel[0,1], self._orbit._satvel[0,2]])
    # 3 x 3 array for starting position storage of actual orient in following order: xx xy xz | yx yy yz | zx zy zz
       self._startV = np.matrix([[1.,0.,0.],
                         [0.,1.,0.],
                         [0.,0.,1.]])
    # 10 x n array for internal storage of actual orient in following order: xx xy xz yx yy yz zx zy zz
       self._orientV =   np.array([[ self._orbit._tins[0],1.,0.,0.,0.,1.,0.,0.,0.,1.],
                                   [ self._orbit._tins[1],1.,0.,0.,0.,1.,0.,0.,0.,1.]])
       self._orientV2 =  np.array([[ self._orbit._tins[0],1.,0.,0.,0.,1.,0.,0.,0.,1.],
                                   [ self._orbit._tins[1],1.,0.,0.,0.,1.,0.,0.,0.,1.]])
       self._orientQ = np.zeros((2, 5))
       self._orientEU= np.zeros((2, 4))
    # 7 x n array for internal storage of actual rotational in deg per second in following order: x y z
       self._rotV =  np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                               [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
       self._rotV2 = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                               [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
    # 7 x n array for internal storage of actual rotational acceleration in deg per square second in following order: x y z
       self._rotaccV  = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                  [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
       self._rotaccV2 = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                  [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
    # 7 x n array for internal storage of actual rotational acceleration rate in deg per square second  and torque momentum in Newton metre | in following order: x y z
       self._torqmomV  = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                   [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
       self._torqmomV2 = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                   [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
    # 7 x n array for internal storage of actual rotational  rate in deg per second  and angular momentum in kilogram square metre per second | in following order: x y z
       self._angmomV =  np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                  [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
       self._angmomV2 = np.array([[ self._orbit._tins[0],0.,0.,0.,0.,0.,0.],
                                  [ self._orbit._tins[1],0.,0.,0.,0.,0.,0.]])
    # implemented as matrix for multipication options
       self._temp1V = np.matrix([[1.,0.,0.],
                         [0.,1.,0.],
                         [0.,0.,1.]])
       self._temp2V = np.matrix([[1.,0.,0.],
                         [0.,1.,0.],
                         [0.,0.,1.]])
       self._anglefreetumble = np.zeros(3)
       self._anglefixspin    = 0.
       #print('Sat Object initialised')
    
    # specifies Orbit by an Orbit Object
    def defOrbit(self, orbloc):
        self._orbit = orbloc
        self._count = self._orbit._count
        # redefines initial rotation matrix in local coords and local variables (global sys as zenit system so - nadir)
        locstartVzenit =  np.array([self._orbit._satv[0,0],self._orbit._satv[0,1],self._orbit._satv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0],self._orbit._satvel[0,1],self._orbit._satvel[0,2]])
        orient = self._q.buildDcmByV(locstartVzenit, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV =  tempstartV * self._rotationM.transpose()
        print('Orbit for Sat defined')        

    # specifies MoI by a MoI Object        
    def defJ(self, jloc):      
        self._J  =   np.matrix([[jloc[0,0],jloc[0,1],jloc[0,2]],
                                [jloc[1,0],jloc[1,1],jloc[1,2]],
                                [jloc[2,0],jloc[2,1],jloc[2,2]]])
        princ, evm = self._jtens.dirCalcEig(jloc)
        self._Jprinc = np.matrix([[princ[0,0],0,0],
                                  [0,princ[1,1],0],
                                  [0,0,princ[2,2]]])
        self._Jevm  =   np.matrix([[evm[0,0],evm[0,1],evm[0,2]],
                                   [evm[1,0],evm[1,1],evm[1,2]],
                                   [evm[2,0],evm[2,1],evm[2,2]]])
        print('Moment of Inertia for Sat defined')
     
    # returns actual Orbit Object
    def getOrb(self):
        return self._orbit
    
    # returns actual MoI Object
    def getJ(self):
        return self._J 
    
    # specifies Nadirs Vectors with Nadir Axe and Flight Axis
    def defNadirV(self, nadirV, flightV):
        # calculate rotational matrix from input (global sys as zenit system so - nadir)
        rot = self._q.buildDcmByV(-nadirV, flightV)
        self._rotationM = np.matrix([[rot[0,0],rot[0,1],rot[0,2]],
                                     [rot[1,0],rot[1,1],rot[1,2]],
                                     [rot[2,0],rot[2,1],rot[2,2]]])
        # redefines initial normalised rotation matrix in local coords and global variables (global sys as zenit system so - nadir)
        locstartVzenit =   np.array([self._orbit._satv[0,0],   self._orbit._satv[0,1],   self._orbit._satv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0], self._orbit._satvel[0,1], self._orbit._satvel[0,2]])
        orient = self._q.buildDcmByV(locstartVzenit, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV =  tempstartV * self._rotationM.transpose()
        print('Starting Vector and Orient Vector redefined for Simple Nadir Pointing Mode')

        
    # simple nadir mode without spin
    def simplenadir(self):
        # reads orbit object for data size and implements array for orientV
        self._orientV = np.zeros((self._count,10))
        self._rotV = np.zeros((self._count-1,7))
        self._temp2V[:,:] = self._startV[:,:]
        for i in range(0, self._count):
            # in local variable
            loctemp1Vnadir  =  np.array([self._orbit._satv[i,0],   self._orbit._satv[i,1],   self._orbit._satv[i,2]])
            loctemp1Vflight =  np.array([self._orbit._satvel[i,0], self._orbit._satvel[i,1], self._orbit._satvel[i,2]])
            orient          =  self._q.buildDcmByV(loctemp1Vnadir, loctemp1Vflight)
            # local var for shorter notation
            temp1orientV    = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                         [orient[1,0],orient[1,1],orient[1,2]],
                                         [orient[2,0],orient[2,1],orient[2,2]]])
            # rotates satellite in orientation in global coords 
            self._temp1V = temp1orientV * self._rotationM.transpose() 
            self._orientV[i,:] = np.array([self._orbit._tins[i], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
            if i > 0:
                # rotational matrix from temp1 to temp2, further calc in local vars
                angleRad, anglelocRad = self._q.getAngleByDcm(self._temp2V,self._temp1V)
                angle    = angleRad    / np.pi * 180 / self._orbit._tincr
                angleloc = anglelocRad / np.pi * 180 / self._orbit._tincr
                self._rotV[i-1,:] = np.array([self._orbit._tins[i-1],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
            self._temp2V[:,:] = self._temp1V[:,:] 
        print('Orientation Vector and Velocity Vector redefined for Simple Nadir Pointing  Mode')
      
        
    # free tumble mode with defined rotation rate about satellite axis (in local coordinates and degree per second)
    def defFreetumbleV(self, nadirV, flightV, alpha, beta, gamma):
        # calculate rotational matrix from input (global sys as zenit system so - nadir)
        rot = self._q.buildDcmByV(-nadirV, flightV)
        self._rotationM   = np.matrix([[rot[0,0],rot[0,1],rot[0,2]],
                                       [rot[1,0],rot[1,1],rot[1,2]],
                                       [rot[2,0],rot[2,1],rot[2,2]]])
        # defines rotation matrix for tumble mode in global coords and for time  in local sat coordinates
        self._anglefreetumble[0] = alpha * self._orbit._tincr / 180 * np.pi
        self._anglefreetumble[1] = beta  * self._orbit._tincr / 180 * np.pi
        self._anglefreetumble[2] = gamma * self._orbit._tincr / 180 * np.pi
        # redefines initial normalised rotation matrix in local coords and global variables (global sys as zenit system so - nadir)
        locstartVzenit  =  np.array([self._orbit._satv[0,0],   self._orbit._satv[0,1],   self._orbit._satv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0], self._orbit._satvel[0,1], self._orbit._satvel[0,2]])
        orient          =  self._q.buildDcmByV(locstartVzenit, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV =  tempstartV * self._rotationM.transpose()
        print('Starting Vector and Orient Vector redefined for Free Tumble Mode')


        
    # free tumble with constant spin rate 
    def freetumble(self):
        # reads orbit object for data size and implements array for orientV
        self._orientV = np.zeros((self._count,10))
        self._rotV  = np.zeros((self._count-1,7))
        # first orientation is similar to start orientation
        self._temp1V[:,:] = self._startV[:,:]
        self._orientV[0,:] = np.array([self._orbit._tins[0], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
        for i in range(0, self._count-1):
            self._temp2V[:,:]  = self._temp1V[:,:] 
            loctemp1V = np.array([[self._temp1V[0,0],self._temp1V[0,1],self._temp1V[0,2]],
                                  [self._temp1V[1,0],self._temp1V[1,1],self._temp1V[1,2]],
                                  [self._temp1V[2,0],self._temp1V[2,1],self._temp1V[2,2]]]) 
            # Rotation is performed in Quaternion class
            loctemp2v, angleRad, anglelocRad = self._q.SatRotFree(self._anglefreetumble, loctemp1V)
            angle    = angleRad    / np.pi * 180 / self._orbit._tincr
            angleloc = anglelocRad / np.pi * 180 / self._orbit._tincr
            # orthogonalise to avoid numerical errors
            loc3  =    np.array([loctemp2v[0,2],loctemp2v[1,2],loctemp2v[2,2]])
            loc1      =    np.array([loctemp2v[0,0],loctemp2v[1,0],loctemp2v[2,0]])
            orient      =    self._q.buildDcmByV(loc3, loc1)
            # local var for shorter notation
            self._temp1V   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                        [orient[1,0],orient[1,1],orient[1,2]],
                                        [orient[2,0],orient[2,1],orient[2,2]]])
            # Writes DCM into Storrage Array
            self._orientV[i+1,:] = np.array([self._orbit._tins[i+1], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
            # Writes local and global rotation rate
            self._rotV[i,:] = np.array([self._orbit._tins[i],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2],])
        print('Orientation Vector and Velocity Vector redefined for Free Tumble Mode')
        
    
    
    # spin nadir mode with defined rotation rate about satellite nadir axis (in local coordinates and degree per second)
    def defNadirspinV(self, nadirV, flightV, gamma):
        # calculate rotational matrix from input (global sys as zenit system so - nadir)
        rot = self._q.buildDcmByV(-nadirV, flightV)
        self._rotationM  = np.matrix([[rot[0,0],rot[0,1],rot[0,2]],
                                      [rot[1,0],rot[1,1],rot[1,2]],
                                      [rot[2,0],rot[2,1],rot[2,2]]])
        # defines rotation matrix for tumble mode in global coords and for time  in local sat coordinates
        self._anglefixspin    = gamma * self._orbit._tincr / 180 * np.pi
        # redefines initial normalised rotation matrix in local coords and global variables
        locstartVzenit  =  np.array([self._orbit._satv[0,0],   self._orbit._satv[0,1],   self._orbit._satv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0], self._orbit._satvel[0,1], self._orbit._satvel[0,2]])
        orient          =  self._q.buildDcmByV(locstartVzenit, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV = tempstartV * self._rotationM.transpose()
        print('Starting Vector and Orient Vector redefined for Nadir Pointing  Spin Mode')
        
        
    
    
    # nadir mode with spin about the nadir axis
    def nadirspin(self):
        # reads orbit object for data size and implements array for orientV
        self._orientV = np.zeros((self._count,10))
        self._rotV = np.zeros((self._count-1,7))
        self._temp1V[:,:] = self._startV[:,:]
        self._orientV[0,:] = np.array([self._orbit._tins[0], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
        for i in range(0, self._count-1):
            self._temp2V[:,:]  = self._temp1V[:,:] 
            loctemp1V = np.array([[self._temp1V[0,0],self._temp1V[0,1],self._temp1V[0,2]],
                                  [self._temp1V[1,0],self._temp1V[1,1],self._temp1V[1,2]],
                                  [self._temp1V[2,0],self._temp1V[2,1],self._temp1V[2,2]]]) 
            # in local variable nadir vectors for continius nadir pointing
            loctemp1Vnadir2  =  np.array([self._orbit._satv[i+1,0],self._orbit._satv[i+1,1],self._orbit._satv[i+1,2]])
            loctemp1Vnadir   =  loctemp1Vnadir2  / np.linalg.norm(loctemp1Vnadir2)
            # calculation in quateul
            self._temp1V, angleRad, anglelocRad =   self._q.SatRotFix(self._anglefixspin, loctemp1V,  self._rotationM, loctemp1Vnadir)
            self._orientV[i+1,:] = np.array([self._orbit._tins[i+1], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
            # rotational matrix from temp1 to temp2, further calc in local vars
            angle    = angleRad    / np.pi * 180 / self._orbit._tincr
            angleloc = anglelocRad / np.pi * 180 / self._orbit._tincr
            self._rotV[i,:] = np.array([self._orbit._tins[i],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
        print('Orientation Vector and Velocity Vector redefined for Nadir Pointing  Spin Mode') 
        
        
    # spin with fixed sun axis mode with defined rotation rate about satellite sun axis (in local coordinates and degree per second)
    def defSunspinV(self, sunV, flightV, gamma):
        # calculate rotational matrix from input (global sys as zenit system so - sunV)
        rot = self._q.buildDcmByV(-sunV, flightV)
        self._rotationM  = np.matrix([[rot[0,0],rot[0,1],rot[0,2]],
                                      [rot[1,0],rot[1,1],rot[1,2]],
                                      [rot[2,0],rot[2,1],rot[2,2]]])
        # defines rotation matrix for tumble mode in global coords and for time  in local sat coordinates
        self._anglefixspin     = gamma * self._orbit._tincr / 180 * np.pi
        # redefines initial normalised rotation matrix in local coords and global variables
        locstartVsun  =    np.array([self._orbit._sunsatv[0,0], self._orbit._sunsatv[0,1],  self._orbit._sunsatv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0],  self._orbit._satvel[0,1],   self._orbit._satvel[0,2]])
        orient          =  self._q.buildDcmByV(locstartVsun, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV = tempstartV * self._rotationM.transpose()
        print('Starting Vector and Orient Vector redefined for Sun Pointing Spin Mode')
        
        
    # sun pointing mode with spin about the nadir axis
    def sunspin(self):
        # reads orbit object for data size and implements array for orientV
        self._orientV = np.zeros((self._count,10))
        self._rotV = np.zeros((self._count-1,7))
        self._temp1V[:,:] = self._startV[:,:]
        self._orientV[0,:] = np.array([self._orbit._tins[0], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
        for i in range(0, self._count-1):
            self._temp2V[:,:]  = self._temp1V[:,:] 
            loctemp1V = np.array([[self._temp1V[0,0],self._temp1V[0,1],self._temp1V[0,2]],
                                  [self._temp1V[1,0],self._temp1V[1,1],self._temp1V[1,2]],
                                  [self._temp1V[2,0],self._temp1V[2,1],self._temp1V[2,2]]]) 
            # in local variable nadir vectors for continius nadir pointing
            loctemp1Vsun2  =  np.array([self._orbit._sunsatv[i+1,0],self._orbit._sunsatv[i+1,1],self._orbit._sunsatv[i+1,2]])
            loctemp1Vsun   = loctemp1Vsun2 / np.linalg.norm(loctemp1Vsun2)
            # calculation in quateul
            self._temp1V, angleRad, anglelocRad =   self._q.SatRotFix(self._anglefixspin, loctemp1V,  self._rotationM, loctemp1Vsun)
            self._orientV[i+1,:] = np.array([self._orbit._tins[i+1], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
            # rotational matrix from temp1 to temp2, further calc in local vars
            angle    = angleRad    / np.pi * 180 / self._orbit._tincr
            angleloc = anglelocRad / np.pi * 180 / self._orbit._tincr
            self._rotV[i,:] = np.array([self._orbit._tins[i],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
        print('Orientation Vector and Velocity Vector redefined for Sun Pointing Spin Mode')
        
    # correction of orientation and spin rate by specified Moment of Inertia (Principal) 
    def SatCorrJ(self):
        # from class MoI with defined _J
        self._Jprinc, self._Jevm = self._jtens.giveMoI()
        self._rotV2 = np.zeros((self._count-1,7))
        self._orientV2 = np.zeros((self._count,10))
        for i in range(0, self._count):
            loctemp2V= np.array([[self._orientV[i,1],self._orientV[i,4],self._orientV[i,7]],
                                 [self._orientV[i,2],self._orientV[i,5],self._orientV[i,8]],
                                 [self._orientV[i,3],self._orientV[i,6],self._orientV[i,9]]])
            # Rotation in Quateul Class
            loctemp1V = self._q.revDcmRotOfDcm(self._Jevm, loctemp2V)
            self._orientV2[i,:] = np.array([self._orbit._tins[i], loctemp1V[0,0], loctemp1V[1,0], loctemp1V[2,0], loctemp1V[0,1], loctemp1V[1,1], loctemp1V[2,1], loctemp1V[0,2], loctemp1V[1,2], loctemp1V[2,2]])
            if i < (self._count-1):
                angle2    = np.array([self._rotV[i,1],self._rotV[i,2],self._rotV[i,3]]) / 180 * np.pi
                angleloc2 = np.array([self._rotV[i,4],self._rotV[i,5],self._rotV[i,6]]) / 180 * np.pi
                # Rotation in Quateul Class
                angle    = self._q.revDcmRotOfEuler(self._Jevm, angle2) / np.pi * 180
                angleloc = self._q.revDcmRotOfEuler(self._Jevm, angleloc2) / np.pi * 180
                self._rotV2[i,:] = np.array([self._orbit._tins[i],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
        print('Satellite Vector and Rotation corrected by Principals of Moment of Inertia')
        
        
    # correction of orientation and spin rate by specified DCM
    def SatCorr(self, dcmRot):
        self._rotV2 = np.zeros((self._count-1,7))
        self._orientV2 = np.zeros((self._count,10))
        for i in range(0, self._count):
            loctemp2V= np.array([[self._orientV[i,1],self._orientV[i,4],self._orientV[i,7]],
                                 [self._orientV[i,2],self._orientV[i,5],self._orientV[i,8]],
                                 [self._orientV[i,3],self._orientV[i,6],self._orientV[i,9]]])
            # Rotation in Quateul Class
            loctemp1V = self._q.revDcmRotOfDcm(dcmRot, loctemp2V)
            self._orientV2[i,:] = np.array([self._orbit._tins[i], loctemp1V[0,0], loctemp1V[1,0], loctemp1V[2,0], loctemp1V[0,1], loctemp1V[1,1], loctemp1V[2,1], loctemp1V[0,2], loctemp1V[1,2], loctemp1V[2,2]])
            if i < (self._count-1):
                angle2    = np.array([self._rotV[i,1],self._rotV[i,2],self._rotV[i,3]]) / 180 * np.pi
                angleloc2 = np.array([self._rotV[i,4],self._rotV[i,5],self._rotV[i,6]]) / 180 * np.pi
                # Rotation in Quateul Class
                angle    = self._q.revDcmRotOfEuler(dcmRot, angle2) / np.pi * 180
                angleloc = self._q.revDcmRotOfEuler(dcmRot, angleloc2) / np.pi * 180
                self._rotV2[i,:] = np.array([self._orbit._tins[i],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
        print('Satellite Vector and Rotation corrected by defined DCM')
        
    
    # acceleration rate for rotation in initial system
    def CalcAccel(self):
        rshape = self._rotV.shape
        count = rshape[0]
        self._rotaccV = np.zeros((count - 1, 7))
        for i in range(0, count - 1):
            anglet1 = np.array([self._rotV[i,1],   self._rotV[i,2],   self._rotV[i,3],   self._rotV[i,4],self._rotV[i,5],     self._rotV[i,6]])
            anglet2 = np.array([self._rotV[i+1,1], self._rotV[i+1,2], self._rotV[i+1,3], self._rotV[i+1,4],self._rotV[i+1,5], self._rotV[i+1,6]])
            angleacc = (anglet2 - anglet1) / self._orbit._tincr 
            self._rotaccV[i,:] = np.array([self._orbit._tins[i],angleacc[0],angleacc[1],angleacc[2],angleacc[3],angleacc[4],angleacc[5]])
        print('Acceleration Rate for defined rotation rate in initial system')
        
    # acceleration rate for rotation in corrected system
    def CalcAccelCorr(self):
        rshape = self._rotV2.shape
        count = rshape[0]
        self._rotaccV2 = np.zeros((count - 1, 7))
        for i in range(0, count - 1):
            anglet1 = np.array([self._rotV2[i,1],   self._rotV2[i,2],   self._rotV2[i,3],   self._rotV2[i,4],   self._rotV2[i,5],  self._rotV2[i,6]])
            anglet2 = np.array([self._rotV2[i+1,1], self._rotV2[i+1,2], self._rotV2[i+1,3], self._rotV2[i+1,4], self._rotV[i+1,5], self._rotV2[i+1,6]])
            angleacc = (anglet2 - anglet1) / self._orbit._tincr 
            self._rotaccV2[i,:] = np.array([self._orbit._tins[i],angleacc[0],angleacc[1],angleacc[2],angleacc[3],angleacc[4],angleacc[5]])
        print('Acceleration Rate for defined rotation rate in initial system')
        
    # angular momentum and torque momentum for initial system    
    def CalcMandL(self):
        rshape = self._rotV.shape
        count = rshape[0]
        rshape2 = self._rotaccV.shape
        count2 = rshape2[0]
        if count2 != count-1:
            self.getAccel()
        J = np.zeros((3,3))
        J[:,:] = self._J[:,:]
        self._torqmomV = np.zeros((count-1, 7))
        self._angmomV = np.zeros((count, 7))
        for i in range(0, count):
            # spin rate for local shortent notation and transformed in rad per second
            omega =  np.array([self._rotV[i,4], self._rotV[i,5],self._rotV[i,6]]) 
            # composing angular momentun
            a0 = (J[0,0] * omega[0] + J[0,1] * omega[1] + J[0,2] * omega[2]) / 180 * np.pi
            a1 = (J[1,0] * omega[0] + J[1,1] * omega[1] + J[1,2] * omega[2]) / 180 * np.pi
            a2 = (J[2,0] * omega[0] + J[2,1] * omega[1] + J[2,2] * omega[2]) / 180 * np.pi
            self._angmomV[i,:] = np.array([self._orbit._tins[i], omega[0], omega[1], omega[2], a0, a1, a2])
            if i < (count-1):
                # spin rate for local shortent notation and transformed in rad per second
                alpha =  np.array([self._rotaccV[i,4], self._rotaccV[i,5],self._rotaccV[i,6]]) 
                # composing angular momentun
                t0 = (J[0,0] * alpha[0] + J[0,1] * alpha[1] + J[0,2] * alpha[2]) / 180 * np.pi
                t1 = (J[1,0] * alpha[0] + J[1,1] * alpha[1] + J[1,2] * alpha[2]) / 180 * np.pi
                t2 = (J[2,0] * alpha[0] + J[2,1] * alpha[1] + J[2,2] * alpha[2]) / 180 * np.pi
                self._torqmomV[i,:] = np.array([self._orbit._tins[i], alpha[0], alpha[1], alpha[2], t0, t1, t2])
        print('Angular Momentum and Torque Momentum  for defined rotation rate in initial system')
        
        
    # angular momentum and torque momentum for corrected system    
    def CalcMandLCorr(self):
        rshape = self._rotV2.shape
        count = rshape[0]
        rshape2 = self._rotaccV2.shape
        count2 = rshape2[0]
        if count2 != count-1:
            self.getAccelCorr()
        J = np.zeros((3,3))
        J[:,:] = self._J[:,:]
        self._torqmomV2 = np.zeros((count-1, 7))
        self._angmomV2 = np.zeros((count, 7))
        for i in range(0, count):
            # spin rate for local shortent notation and transformed in rad per second
            omega =  np.array([self._rotV2[i,4], self._rotV2[i,5],self._rotV2[i,6]]) 
            # composing angular momentun
            a0 = (J[0,0] * omega[0] + J[0,1] * omega[1] + J[0,2] * omega[2]) / 180 * np.pi
            a1 = (J[1,0] * omega[0] + J[1,1] * omega[1] + J[1,2] * omega[2]) / 180 * np.pi
            a2 = (J[2,0] * omega[0] + J[2,1] * omega[1] + J[2,2] * omega[2]) / 180 * np.pi
            self._angmomV2[i,:] = np.array([self._orbit._tins[i], omega[0], omega[1], omega[2], a0, a1, a2])
            if i < (count-1):
                # spin rate for local shortent notation and transformed in rad per second
                alpha =  np.array([self._rotaccV2[i,4], self._rotaccV2[i,5],self._rotaccV2[i,6]]) 
                # composing angular momentun
                t0 = (J[0,0] * alpha[0] + J[0,1] * alpha[1] + J[0,2] * alpha[2]) / 180 * np.pi
                t1 = (J[1,0] * alpha[0] + J[1,1] * alpha[1] + J[1,2] * alpha[2]) / 180 * np.pi
                t2 = (J[2,0] * alpha[0] + J[2,1] * alpha[1] + J[2,2] * alpha[2]) / 180 * np.pi
                self._torqmomV2[i,:] = np.array([self._orbit._tins[i], alpha[0], alpha[1], alpha[2], t0, t1, t2])
        print('Angular Momentum and Torque Momentum  for defined rotation rate in initial system')
        
        
    # returns satellite orientation and vector in DCM
    def getSatDCM(self):
        return self._orientV
    
    # returns satellite orientation and vector in Quaternion
    def getSatQuaternion(self):
        self._orientQ = np.zeros((self._count,5))
        for i in range(0, self._count):
            tempdcm = np.matrix([[self._orientV[i,1],self._orientV[i,4],self._orientV[i,7]],
                                 [self._orientV[i,2],self._orientV[i,5],self._orientV[i,8]],
                                 [self._orientV[i,3],self._orientV[i,6],self._orientV[i,9]]])
            tempquat = self._q.dirDcmToQuat(tempdcm)
            self._orientQ[i,:] = np.array([self._orientV[i,0],tempquat[0],tempquat[1],tempquat[2],tempquat[3]])
        return self._orientQ
    
    # returns satellite orientation and vector in Euler
    def getSatEuler(self):
        self._orientEU = np.zeros((self._count,4))
        for i in range(0, self._count):
            tempdcm = np.matrix([[self._orientV[i,1],self._orientV[i,4],self._orientV[i,7]],
                                 [self._orientV[i,2],self._orientV[i,5],self._orientV[i,8]],
                                 [self._orientV[i,3],self._orientV[i,6],self._orientV[i,9]]])
            tempeuler = self._q.dirDcmToEuler(tempdcm)
            self._orientEU[i,:] = np.array([self._orientV[i,0],tempeuler[0],tempeuler[1],tempeuler[2]])
        return self._orientEU
    
    
    
    # specifies Point Vectors with Point Axe and Flight Axis
    def defNadirPointV(self, nadirV, flightV):
        # calculate rotational matrix from input (global sys as zenit system so - nadir)
        rot = self._q.buildDcmByV(-nadirV, flightV)
        self._rotationM = np.matrix([[rot[0,0],rot[0,1],rot[0,2]],
                                     [rot[1,0],rot[1,1],rot[1,2]],
                                     [rot[2,0],rot[2,1],rot[2,2]]])
        # redefines initial normalised rotation matrix in local coords and global variables (global sys as zenit system so - nadir)
        locstartVzenit =   np.array([self._orbit._earthLocPosSatv[0,0],   self._orbit._earthLocPosSatv[0,1],   self._orbit._earthLocPosSatv[0,2]])
        locstartVflight =  np.array([self._orbit._satvel[0,0], self._orbit._satvel[0,1], self._orbit._satvel[0,2]])
        orient = self._q.buildDcmByV(locstartVzenit, locstartVflight)
        # local var for shorter notation
        tempstartV   = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                  [orient[1,0],orient[1,1],orient[1,2]],
                                  [orient[2,0],orient[2,1],orient[2,2]]])
        # rotates satellite in initial orientation in global coords 
        self._startV =  tempstartV * self._rotationM.transpose()
        print('Starting Vector and Orient Vector redefined for Nadir Pointing Mode')

        
    #  nadir pointing mode without spin
    def nadirpoint(self):
        # reads orbit object for data size and implements array for orientV
        self._orientV = np.zeros((self._count,10))
        self._rotV = np.zeros((self._count-1,7))
        self._temp2V[:,:] = self._startV[:,:]
        for i in range(0, self._count):
            # in local variable
            loctemp1Vnadir  =  np.array([self._orbit._earthLocPosSatv[i,0],   self._orbit._earthLocPosSatv[i,1],   self._orbit._earthLocPosSatv[i,2]])
            loctemp1Vflight =  np.array([self._orbit._satvel[i,0], self._orbit._satvel[i,1], self._orbit._satvel[i,2]])
            orient          =  self._q.buildDcmByV(loctemp1Vnadir, loctemp1Vflight)
            # local var for shorter notation
            temp1orientV    = np.matrix([[orient[0,0],orient[0,1],orient[0,2]],
                                         [orient[1,0],orient[1,1],orient[1,2]],
                                         [orient[2,0],orient[2,1],orient[2,2]]])
            # rotates satellite in orientation in global coords 
            self._temp1V = temp1orientV * self._rotationM.transpose() 
            self._orientV[i,:] = np.array([self._orbit._tins[i], self._temp1V[0,0],self._temp1V[1,0],self._temp1V[2,0],self._temp1V[0,1],self._temp1V[1,1],self._temp1V[2,1],self._temp1V[0,2],self._temp1V[1,2],self._temp1V[2,2]])
            if i > 0:
                # rotational matrix from temp1 to temp2, further calc in local vars
                angleRad, anglelocRad = self._q.getAngleByDcm(self._temp2V,self._temp1V)
                angle    = angleRad    / np.pi * 180 / self._orbit._tincr
                angleloc = anglelocRad / np.pi * 180 / self._orbit._tincr
                self._rotV[i-1,:] = np.array([self._orbit._tins[i-1],angle[0],angle[1],angle[2],angleloc[0],angleloc[1],angleloc[2]])
            self._temp2V[:,:] = self._temp1V[:,:] 
        print('Orientation Vector and Velocity Vector redefined for Nadir Pointing  Mode')
        
  
      
 