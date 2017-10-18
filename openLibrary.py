# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 13:35:17 2017
Libary Class
@author: Danilo
"""
import numpy as np
from geosat import SurfaceGeometry as sfg
from dissipation import *
from GeometricShader import *

    #simple CubeSat 1 U all surfaces with 
class Cube:
    
    def __init__(self):
        self._count = 6
        # no array but a list
        self._sf = [sfg(),sfg(),sfg(),sfg(),sfg(),sfg()]
        self._orient = np.array([[ 1,  0,  0],
                        [ 0,  1,  0],
                        [-1,  0,  0],
                        [ 0, -1,  0],
                        [ 0,  0,  1],
                        [ 0,  0, -1]])
        self._name = 'surface'
        surface = SolarBoard1U()
        for i in range(0, self._count):
            surface.setSurface(self._sf[i])
            self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        print('Simple Cube Sat Initialised')
        
        # payload optic without deployed solar 2U arrays
class Cube2U:

    
    def __init__(self):
        self._count = 7
        # no array but a list
        self._sf = [sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg()]
        self._orient = np.array([[ 1,  0,  0],
                        [ 0,  1,  0],
                        [-1,  0,  0],
                        [-1,  0,  0],
                        [ 0, -1,  0],
                        [ 0,  0,  1],
                        [ 0,  0, -1]])
        self._name = 'surface'
        surface1U = SolarBoard1U()
        surface2U = SolarBoard2U()
        surface1UspecSolarcell = SolarBoard1USpecial()
        surface1UspecPayload   = Board1UPayload()
        i = 0
        surface2U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 1
        surface2U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 2
        surface1UspecSolarcell.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 3
        surface1UspecPayload.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 4
        surface2U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 5
        surface1U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 6
        surface1U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        print('2U Cube Sat Initialised')
        
        
                # payload optic with deployed solar 2U arrays
class Cube2Udeploy:

    
    def __init__(self):
        self._count = 9
        # no array but a list
        self._sf = [sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg()]
        self._orient = np.array([[ 1,  0,  0],
                                 [-1,  0,  0],
                                 [ 0, -1,  0],
                                 [-1,  0,  0],
                                 [-1,  0,  0],
                                 [ 0,  1,  0],
                                 [-1,  0,  0],
                                 [ 0,  0,  1],
                                 [ 0,  0, -1]])
    # down solar polar -x
        self._dcm1  = np.array([[ 0,0,-1],
                                [-1,0, 0],
                                [ 0,1, 0]])

    # down cube -y
        self._dcm2  = np.array([[-1, 0, 0],
                                [ 0, 0,-1],
                                [ 0,-1, 0]])



    # down cube y
        self._dcm3  = np.array([[-1,0,0],
                                [0,0,1],
                                [0,1,0]])
    # down solar polar -x
        self._dcm4  = np.array([[ 0, 0,-1],
                                [ 1, 0, 0],
                                [ 0,-1, 0]])
    
        self._name = 'surface'
    # shader
        self._length = 100
        self._height = 100
        self._pos    = 0
        self._areawidth = 200
        self._objectwidth = 200
        surface1U = SolarBoard1U()
        surface2U = SolarBoard2U()
        surface6U = SolarBoard6U()
        surface2Ualu = Alu2U()
        surface1UspecSolarcell = SolarBoard1USpecial()
        surface1UspecPayload   = Board1UPayload()
        i = 0
        surface6U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 1
        surface2Ualu.setSurface(self._sf[i])
        shader1 = ObjectShader(self._length, self._height, self._areawidth, self._objectwidth, self._pos)
        self._sf[i].initShadedSurface(self._dcm1, self._name, (i + 1), shader1)
        i = 2
        surface2U.setSurface(self._sf[i])
        shader2 = ObjectShader(self._length, self._height, self._areawidth, self._objectwidth, self._pos)
        self._sf[i].initShadedSurface(self._dcm2, self._name, (i + 1), shader2)
        i = 3
        surface1UspecSolarcell.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 4
        surface1UspecPayload.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 5
        surface2U.setSurface(self._sf[i])
        shader3 = ObjectShader(self._length, self._height, self._areawidth, self._objectwidth, self._pos)
        self._sf[i].initShadedSurface(self._dcm3, self._name, (i + 1), shader3)
        i = 6
        surface2Ualu.setSurface(self._sf[i])
        shader4 = ObjectShader(self._length, self._height, self._areawidth, self._objectwidth, self._pos)
        self._sf[i].initShadedSurface(self._dcm4, self._name, (i + 1), shader4)
        i = 7
        surface1U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        i = 8
        surface1U.setSurface(self._sf[i])
        self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        print('2U Cube Sat with deployed Solararrays Initialised')
        
        
        
class Octa:
    
    def __init__(self):
        self._count = 10
    # no array but a list
        self._sf = [sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg(),sfg()]
        self._orient = np.array([[ 1,  0,  0],
                        [ 1,  1,  0],
                        [ 0,  1,  0],
                        [-1,  1,  0],
                        [-1,  0,  0],
                        [-1, -1,  0],
                        [ 0, -1,  0],
                        [ 1, -1,  0],
                        [ 0,  0,  1],
                        [ 0,  0, -1]])
        self._name = 'surface'
        for i in range(0, self._count):
            self._sf[i].initSimpleSurface(self._orient[i,:], self._name, (i + 1))
        print('Simple Cube Sat Initialised')
        
        
class SolarCell1:
    # Standard Cell Size Auzurspace 3G30C Advanced
    _areasize   = 3000
    _efficiency = 0.25
    _alpha      = 0.91
    _epsilon    = 0.78
    
class SolarCell2:
    # Standard Cell Size Auzurspace 3G30C Advanced
    _areasize   = 375
    _efficiency = 0.25
    _alpha      = 0.91
    _epsilon    = 0.78
    
    
    # Optical Parameters
class SurfaceOpen:
   _alphaPCB    = 0.5
   _epsilonPCB  = 0.9
   _alphaAlu    = 0.192
   _epsilonAlu  = 0.245
   _alphaDNC   = 0.453
   _epsilonDNC = 0.132
   _alphaPatch   = 0.346
   _epsilonPatch = 0.028

    
    
  # solar Board for CubeSat     
class SolarBoard1U:


        # init Solarboard 
    def __init__(self):
            # Standard Solarboard 1u Cube
        self._areasize    = 100 * 100
        self._sizePCBges  = 95  *  82
    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges
        self._Solararea   = SolarCell1._areasize * 2
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaAlu    +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonAlu  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell1._alpha, SolarCell1._epsilon, SolarCell1._efficiency)
 
  # special solar Board for CubeSat     
class SolarBoard1USpecial:
    

        # init Solarboard 
    def __init__(self):
    # Special Solarboard 1U
        self._areasize    = 100 * 100
        self._sizePCBges  = 82  *  80
        self._sizePatch   = 30  *  30
    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges -     self._sizePatch
        self._Solararea   = SolarCell2._areasize * 8
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._sizePatch * SurfaceOpen._alphaPatch   + self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaDNC    +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._sizePatch * SurfaceOpen._epsilonPatch + self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonDNC  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell2._alpha, SolarCell2._epsilon, SolarCell1._efficiency)

  # special solar Board for CubeSat     
class Board1UPayload:
  

        # init Solarboard 
    def __init__(self):
    # Standard Payloadboard 1U
        self._areasize    = 100 * 100
        self._sizePCBges  = 0

    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges 
        self._Solararea   = SolarCell1._areasize * 0
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaDNC    +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonDNC  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell1._alpha, SolarCell1._epsilon, SolarCell1._efficiency)
     
    
  # solar Board 2U for CubeSat     
class SolarBoard2U:

        # init Solarboard 
    def __init__(self):
            # Standard Solarboard 2U
        self._areasize    = 2 * 100 * 100
        self._sizePCBges  = 2 * 95  *  82
    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges
        self._Solararea   = SolarCell1._areasize * 4
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaDNC   +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonDNC  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell1._alpha, SolarCell1._epsilon, SolarCell1._efficiency)


  # solar Board 2U for CubeSat     
class Alu2U:
 

        # init Solarboard 
    def __init__(self):
    # Standard Aluboard 2U (backside of deployed solararray)
        self._areasize    = 2 * 100 * 100
        self._sizePCBges  = 2 * 95  *  82
    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges
        self._Solararea   = SolarCell1._areasize * 0
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaDNC    +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonDNC  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell1._alpha, SolarCell1._epsilon, SolarCell1._efficiency)


  # solar Board for CubeSat 2 U with deployed solar arrays 
class SolarBoard6U:
    
        # init Solarboard  6U
    def __init__(self):
            # Standard Solarboard 6U
        self._areasize    = 6 * 100 * 100
        self._sizePCBges  = 6 * 95  *  82
    # build Alpha Epsilon
        self._areasizeAlu =     self._areasize -     self._sizePCBges
        self._Solararea   = SolarCell1._areasize * 12
        self._areasizePCB =     self._sizePCBges -     self._Solararea 
        self._alphamix   = ( self._areasizePCB * SurfaceOpen._alphaPCB   +  self._areasizeAlu  * SurfaceOpen._alphaDNC    +  self._Solararea * SolarCell1._alpha   ) / self._areasize
        self._epsilonmix = ( self._areasizePCB * SurfaceOpen._epsilonPCB +  self._areasizeAlu  * SurfaceOpen._epsilonDNC  +  self._Solararea * SolarCell1._epsilon ) / self._areasize

        # redefines handed surface geometry new with sideframe properties
    def setSurface(self, sf):
        sf.setSurfaceParamSolar(self._alphamix, self._epsilonmix, self._areasize, self._Solararea, SolarCell1._alpha, SolarCell1._epsilon, SolarCell1._efficiency)
 

    
    # simple Elecsystem with out solarcells and charge
class SimpleElectroSystem:
    
    # defines all components based on orbit and dissipation
    def __init__(self, orbit):
        self.initComp = initialComponent()
        self.initComp.defineComponent(orbit)
        self.OBC  =         threeModeComponent(self.initComp._count)
        self.PCU  =         simplePCU(self.initComp._count)
        self.COMUHF     =   threeModeComponent(self.initComp._count)
        self.COMSBAND   =   threeModeComponent(self.initComp._count)
        self.AOCS =         twoModeComponent(self.initComp._count)
        self.EBOX =         componentContainer(self.initComp._count, 2)
        self.OBC.defineFactors(0.3, 2., 3., 0, 0)
        self.PCU.defineFactors(0.5, 0.2)
        self.COMUHF.defineFactors(0.0, 0.6, 3.5, 0.35, 0.4)
        self.COMSBAND.defineFactors(0.0, 1.5, 6.5, 0.15, 0.15)
        self.AOCS.defineFactors(0.2, 2.1, 0.3)
        # defines further components, which will hand over to FEM
        self.BAT.compContain.setName('battery')
        self.BAT.compContain.setArea(80*80)
        self.EBOX.setName('ebox')
        self.EBOX.setArea(95*95)
        self.COMSBAND.compContain.setName('sband')
        self.COMSBAND.compContain.setArea(90*45)
        self.sumAll = componentContainer(self.initComp._count, 3) 
        
        # simple example to init a electronic mode controled system without much effort
    def calcDiss(self, timestep, mode):
            # waiting mode
        if mode == 0:
            self.COMUHF.calcDiss(1, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
            # receive mode
        elif mode == 1:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(1, timestep)
            self.OBC.calcDiss(1, timestep)
            self.AOCS.calcDiss(1, timestep)
             # sent mode
        elif mode == 2:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(2, timestep)
            self.OBC.calcDiss(2, timestep)
            self.AOCS.calcDiss(1, timestep)
            # safe mode
        else:
            self.COMUHF.calcDiss(0, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
        # calc PCU depend on other systems
        PCUload = self.AOCS.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.COMSBAND.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        self.PCU.calcDiss(PCUload, timestep)
        # sums all subsystems stored in the ebox
        EboxEconsum  = self.AOCS.compContain.getEconsum(timestep) + self.PCU.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        EboxDiss     = self.AOCS.compContain.getDiss(timestep)    + self.PCU.compContain.getDiss(timestep)    + self.COMUHF.compContain.getDiss(timestep)    + self.OBC.compContain.getDiss(timestep)
        EboxValue    = np.array([EboxEconsum, EboxDiss])
        self.EBOX.writeValue(EboxValue, timestep)
        # sums electric energy consumed in all subsystems in satellite
        sumEconsumAll      =  EboxEconsum + self.COMSBAND.compContain.getEconsum(timestep)
        # heat disstributed in the total system
        sumDissiAll        =  EboxDiss    + self.COMSBAND.compContain.getDiss(timestep) 
        # vector has to be written in rows
        value = np.array([sumEconsumAll, sumDissiAll, mode])
        self.sumAll.writeValue(value, timestep)
        
        
        # Dissipation Power has to be performed after calcEconsum and is needed for further FEM Calc
    def getDissipation(self):
        plotmatrix = np.zeros([self.initComp._count, 2])
        #plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,0] = self.EBOX.container[:,1]
        plotmatrix[:,1] = self.COMSBAND.compContain.container[:,1]
        return plotmatrix
    
            # Heatflux calc from Dissipation and specified Area [mm^2] but returns in [W/m^2] has to be performed after calcEconsum and is needed for further FEM Calc
    def getHeatflux(self):
        plotmatrix = np.zeros([self.initComp._count, 2])
        #plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,0] = self.EBOX.container[:,1] / self.EBOX.getArea() * 10.**6
        plotmatrix[:,1] = self.COMSBAND.compContain.container[:,1] / self.COMSBAND.compContain.getArea() * 10.**6
        return plotmatrix
    
    
        # returns name for further FEM Calc (Output)
    def getName(self):
        name = [self.EBOX.getName(),  self.COMSBAND.compContain.getName()]
        return name
       
          # get plotmatrix for overview
    def getPlot(self):   
        plotmatrix = np.zeros([self.initComp._count, 14])
        plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,1] = self.sumAll.container[:,0]
        plotmatrix[:,2] = self.sumAll.container[:,1]
        plotmatrix[:,3] = self.sumAll.container[:,2]
        plotmatrix[:,4] = self.PCU.compContain.container[:,0]
        plotmatrix[:,5] = self.PCU.compContain.container[:,1]
        plotmatrix[:,6] = self.OBC.compContain.container[:,0]
        plotmatrix[:,7] = self.OBC.compContain.container[:,1]
        plotmatrix[:,8] = self.AOCS.compContain.container[:,0]
        plotmatrix[:,9] = self.AOCS.compContain.container[:,1]
        plotmatrix[:,10] = self.COMUHF.compContain.container[:,0]
        plotmatrix[:,11] = self.COMUHF.compContain.container[:,1]
        plotmatrix[:,12] = self.COMSBAND.compContain.container[:,0]
        plotmatrix[:,13] = self.COMSBAND.compContain.container[:,1]
        return plotmatrix
    
    
    
    
    
    
    # simple Elecsystem with  solarcells and charge
class SimpleChargeElectroSystem:
    
    # defines all components based on scripted dissipation and orbit and surface geometry
    def __init__(self, surfaceheatflux, surfaceGeometry):
        self.initComp = initialComponent()
        self.initComp.defineComponent(surfaceheatflux)
        self.sHeatflux = surfaceheatflux
        # fills list with simpleSolar for better storag, corressponding to the handed surfaceGeometry
        self.SOLAR         = [simpleSolar()]
        self.faceCount     = len(surfaceGeometry) 
        for i in range(0, (self.faceCount - 1)):
            self.SOLAR.append(simpleSolar())    
        self.OBC  =         threeModeComponent(self.initComp._count)
        self.PCU  =         simplePCU(self.initComp._count)
        self.COMUHF     =   threeModeComponent(self.initComp._count)
        self.COMSBAND   =   threeModeComponent(self.initComp._count)
        self.AOCS =         twoModeComponent(self.initComp._count)
        self.EBOX =         componentContainer(self.initComp._count, 2)
        self.BAT  =         simpleBattery(self.initComp._count, self.initComp._tincr)
        self.SUMSOLAR =     componentContainer(self.initComp._count, 3)
        self.OBC.defineFactors(0.3, 2., 3., 0, 0)
        self.COMUHF.defineFactors(0.0, 0.6, 3.5, 0.35, 0.4)
        self.COMSBAND.defineFactors(0.0, 1.5, 6.5, 0.15, 0.15)
        self.PCU.defineFactors(0.5, 0.2)
        self.AOCS.defineFactors(0.2, 2.1, 0.3)
        # defines further components, which will hand over to FEM
        self.BAT.compContain.setName('battery')
        self.BAT.compContain.setArea(80*80)
        self.EBOX.setName('ebox')
        self.EBOX.setArea(95*95)
        self.COMSBAND.compContain.setName('sband')
        self.COMSBAND.compContain.setArea(90*45)
                # defines alla Battery Values (basicDissipation, chargeFactor, dechargeFactor, capacity [Wh], initialChargeFactor [Factor]):
        self.BAT.defineFactors(0.001, 0.1, 0.12, 9.1, 0.9995)
            # redifnes SimpleSolar Generators
        for i in range(0, self.faceCount):
            self.SOLAR[i].defineSolarGenerator(surfaceGeometry[i], self.initComp._count)
        self.sumAll = componentContainer(self.initComp._count, 4) 
        
        
        # simple example to init a electronic mode controled system without much effort
    def calcEconsum(self, timestep, mode):
            # waiting mode
        if mode == 0:
            self.COMUHF.calcDiss(1, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
            # receive mode
        elif mode == 1:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(1, timestep)
            self.OBC.calcDiss(1, timestep)
            self.AOCS.calcDiss(1, timestep)
             # sent mode
        elif mode == 2:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(2, timestep)
            self.OBC.calcDiss(2, timestep)
            self.AOCS.calcDiss(1, timestep)
            # safe mode
        else:
            self.COMUHF.calcDiss(0, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
        # calc PCU depend on other systems
        PCUload = self.AOCS.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.COMSBAND.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        self.PCU.calcDiss(PCUload, timestep)
        # sums all subsystems stored in the ebox
        EboxEconsum  = self.AOCS.compContain.getEconsum(timestep) + self.PCU.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        EboxDiss     = self.AOCS.compContain.getDiss(timestep)    + self.PCU.compContain.getDiss(timestep)    + self.COMUHF.compContain.getDiss(timestep)    + self.OBC.compContain.getDiss(timestep)
        EboxValue    = np.array([EboxEconsum, EboxDiss])
        self.EBOX.writeValue(EboxValue, timestep)
        sumEconsum = EboxEconsum + self.COMSBAND.compContain.getEconsum(timestep) 
        # calculation of energy generation for all surfaces based of incomming heatfluxes
        # vector has to be written in rows
        tempSumSOLAR    = np.zeros(3)
        for i in range(0, self.faceCount):
            tempheatflux   = self.sHeatflux._heatsun[timestep, i] + self.sHeatflux._heatalb[timestep, i] 
            self.SOLAR[i].calcEnergy(tempheatflux, timestep)
            tempSumSOLAR[0] = tempSumSOLAR[0] + self.SOLAR[i].compContain.getMainValue(timestep)
            tempSumSOLAR[2] = tempSumSOLAR[2] + self.SOLAR[i].compContain.getTerValue(timestep)
        # calculates of charge for batteries
        charge = tempSumSOLAR[0] - sumEconsum
        self.BAT.calcDissCharge(charge, timestep)
        # sums all subsystems in satellite
        # this is needed to calculate the amount of the energy transformed into electricity at the solar generators
        sumEconsumAll      =  EboxEconsum + self.BAT.compContain.getEconsum(timestep) + self.COMSBAND.compContain.getEconsum(timestep)
        # only respects energy lost, means energy provided by batteries is not included in calc
        sumEconsumAllSys   =  EboxEconsum + self.BAT.compContain.getDiss(timestep) + self.COMSBAND.compContain.getEconsum(timestep)
        # heat disstributed in the total system
        sumDissiAll        =  EboxDiss    + self.BAT.compContain.getDiss(timestep) + self.COMSBAND.compContain.getDiss(timestep) 
        # vector has to be written in rows
        value = np.array([sumEconsumAll, sumDissiAll, sumEconsumAllSys, mode])
        # calcs soalr energy usage
        tempSumSOLAR[1] = self.SOLAR[0].dirCalcConsumption(sumEconsumAllSys, self.SUMSOLAR.getMainValue())
        solarloadfaktor = tempSumSOLAR[1] / self.SUMSOLAR.getMainValue()
        for i in range(0, self.faceCount):
            self.SOLAR[i].CalcConsumption(solarloadfaktor)
        self.SUMSOLAR.writeValue(tempSumSOLAR, timestep)
        self.sumAll.writeValue(value, timestep)
        
          # get plotmatrix for overview
    def getPlot(self):   
        plotmatrix = np.zeros([self.initComp._count, 20])
        plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,1] = self.sumAll.container[:,0]
        plotmatrix[:,2] = self.sumAll.container[:,1]
        plotmatrix[:,3] = self.sumAll.container[:,2]
        plotmatrix[:,4] = self.sumAll.container[:,3]
        plotmatrix[:,5] = self.BAT.compContain.container[:,0]
        plotmatrix[:,6] = self.BAT.compContain.container[:,1]
        plotmatrix[:,7] = self.SUMSOLAR.container[:,0]
        plotmatrix[:,8] = self.SUMSOLAR.container[:,1]
        plotmatrix[:,9] = self.SUMSOLAR.container[:,2]
        plotmatrix[:,10] = self.PCU.compContain.container[:,0]
        plotmatrix[:,11] = self.PCU.compContain.container[:,1]
        plotmatrix[:,12] = self.OBC.compContain.container[:,0]
        plotmatrix[:,13] = self.OBC.compContain.container[:,1]
        plotmatrix[:,14] = self.AOCS.compContain.container[:,0]
        plotmatrix[:,15] = self.AOCS.compContain.container[:,1]
        plotmatrix[:,16] = self.COMUHF.compContain.container[:,0]
        plotmatrix[:,17] = self.COMUHF.compContain.container[:,1]
        plotmatrix[:,18] = self.COMSBAND.compContain.container[:,0]
        plotmatrix[:,19] = self.COMSBAND.compContain.container[:,1]
        return plotmatrix
    
                  # get solar load for further heat reduction
    def getSolarload(self):
        return self.SUMSOLAR.container[:,1]
    
                # get dissipation for further thermal simuolation
    def getSumDissipation(self):
        return self.sumAll.container[:,1]
    
    
        # Dissipation Power has to be performed after calcEconsum and is needed for further FEM Calc
    def getDissipation(self):
        plotmatrix = np.zeros([self.initComp._count, 3])
        #plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,0] = self.EBOX.container[:,1]
        plotmatrix[:,1] = self.BAT.compContain.container[:,1]
        plotmatrix[:,2] = self.COMSBAND.compContain.container[:,1]
        return plotmatrix
    
            # Heatflux calc from Dissipation and specified Area [mm^2] but returns in [W/m^2] has to be performed after calcEconsum and is needed for further FEM Calc
    def getHeatflux(self):
        plotmatrix = np.zeros([self.initComp._count, 3])
        #plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,0] = self.EBOX.container[:,1] / self.EBOX.getArea() * 10.**6
        plotmatrix[:,1] = self.BAT.compContain.container[:,1] / self.BAT.compContain.getArea() * 10.**6
        plotmatrix[:,2] = self.COMSBAND.compContain.container[:,1] / self.COMSBAND.compContain.getArea() * 10.**6
        return plotmatrix
    
    
        # returns name for further FEM Calc (Output)
    def getName(self):
        name = [self.EBOX.getName(), self.BAT.compContain.getName(), self.COMSBAND.compContain.getName()]
        return name    
    
    
    # Standard Elecsystem with  solarcells and charge
class StandardChargeElectroSystem:
    
    # defines all components based on scripted dissipation and orbit and surface geometry
    def __init__(self, surfaceheatflux):
        self.initComp = initialComponent()
        self.initComp.defineComponent(surfaceheatflux)
        self.sHeatflux = surfaceheatflux
        # fills list with simpleSolar for better storag, corressponding to the handed surfaceGeometry
        self.SOLAR =        standardSolar()
        self.OBC  =         threeModeComponent(self.initComp._count)
        self.PCU  =         simplePCU(self.initComp._count)
        self.COMUHF     =   threeModeComponent(self.initComp._count)
        self.COMSBAND   =   threeModeComponent(self.initComp._count)
        self.AOCS =         twoModeComponent(self.initComp._count)
        self.EPS  =         twoModeComponent(self.initComp._count)
        self.EBOX =         componentContainer(self.initComp._count, 2)
        self.BAT  =         simpleBattery(self.initComp._count, self.initComp._tincr)
        self.OBC.defineFactors(0.1, 0.15, 0.25, 0, 0)
        self.EPS.defineFactors(0.1, 0.25, 0.25)
        self.COMUHF.defineFactors(0.0, 0.6, 2.0, 0.35, 0.4)
        self.COMSBAND.defineFactors(0.0, 1.5, 6.5, 0.15, 0.15)
        self.PCU.defineFactors(0.01, 0.2)
        self.AOCS.defineFactors(0.2, 1.5, 0.3)
                # defines alla Battery Values (basicDissipation, chargeFactor, dechargeFactor, capacity [Wh], initialChargeFactor [Factor]):
        self.BAT.defineFactors(0.001, 0.1, 0.12, 2 * 9.1, 0.98)
            # redifnes SimpleSolar Generators with surfaceheatflux object
            
         # array for compContainers for shortent notation
        self.dissiContain = []
        # defines further components, which will hand over to FEM
        self.BAT.compContain.setName('battery')
        self.BAT.compContain.setArea(80*80)
        self.dissiContain.append(self.BAT.compContain)
        self.EBOX.setName('ebox')
        self.EBOX.setArea(95*95)
        self.dissiContain.append(self.EBOX)
        self.COMSBAND.compContain.setName('sband')
        self.COMSBAND.compContain.setArea(90*45)
        self.dissiContain.append(self.COMSBAND.compContain)
        
        self.count = len(self.dissiContain)
                # defines alla Battery Values (basic
        self.SOLAR.defineSolarGenerator(surfaceheatflux)
        self.sumAll = componentContainer(self.initComp._count, 4) 
        
        
        # simple example to init a electronic mode controled system without much effort
    def calcEconsum(self, timestep, mode):
        # EPS must work
        self.EPS.calcDiss(1, timestep)
            # waiting mode
        if mode == 0:
            self.COMUHF.calcDiss(1, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
            # receive mode
        elif mode == 1:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(1, timestep)
            self.OBC.calcDiss(1, timestep)
            self.AOCS.calcDiss(1, timestep)
             # sent mode
        elif mode == 2:
            self.COMUHF.calcDiss(1, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(1, timestep)
            self.AOCS.calcDiss(1, timestep)
        elif mode == 3:
            self.COMUHF.calcDiss(2, timestep)
            self.COMSBAND.calcDiss(2, timestep)
            self.OBC.calcDiss(2, timestep)
            self.AOCS.calcDiss(1, timestep)
            # safe mode
        else:
            self.COMUHF.calcDiss(0, timestep)
            self.COMSBAND.calcDiss(0, timestep)
            self.OBC.calcDiss(0, timestep)
            self.AOCS.calcDiss(0, timestep)
        # calc PCU depend on other systems
        PCUload = self.EPS.compContain.getEconsum(timestep) + self.AOCS.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.COMSBAND.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        self.PCU.calcDiss(PCUload, timestep)
        # sums all subsystems stored in the ebox
        EboxEconsum  = self.EPS.compContain.getEconsum(timestep) + self.AOCS.compContain.getEconsum(timestep) + self.PCU.compContain.getEconsum(timestep) + self.COMUHF.compContain.getEconsum(timestep) + self.OBC.compContain.getEconsum(timestep)
        EboxDiss     = self.EPS.compContain.getDiss(timestep)    + self.AOCS.compContain.getDiss(timestep)    + self.PCU.compContain.getDiss(timestep)    + self.COMUHF.compContain.getDiss(timestep)    + self.OBC.compContain.getDiss(timestep)
        EboxValue    = np.array([EboxEconsum, EboxDiss])
        self.EBOX.writeValue(EboxValue, timestep)
        sumEconsum = EboxEconsum + self.COMSBAND.compContain.getEconsum(timestep)          
        # calculates of charge for batteries
        charge = self.SOLAR.compContain.getMainValue(timestep) - sumEconsum
        self.BAT.calcDissCharge(charge, timestep)
        # sums all subsystems in satellite
        # this is needed to calculate the amount of the energy transformed into electricity at the solar generators
        sumEconsumAll      =  EboxEconsum + self.BAT.compContain.getEconsum(timestep) + self.COMSBAND.compContain.getEconsum(timestep)
        # only respects energy lost, means energy provided by batteries is not included in calc
        sumEconsumAllSys   =  EboxEconsum + self.BAT.compContain.getDiss(timestep) + self.COMSBAND.compContain.getEconsum(timestep)
        # heat disstributed in the total system
        sumDissiAll        =  EboxDiss    + self.BAT.compContain.getDiss(timestep) + self.COMSBAND.compContain.getDiss(timestep)        
        # vector has to be written in rows
        value = np.array([sumEconsumAll, sumDissiAll, sumEconsumAllSys, mode])
        # calculation of energy generation for all surfaces based of incomming heatfluxes
        self.SOLAR.calcConsumption(sumEconsumAll, timestep)  
        self.sumAll.writeValue(value, timestep)
        
          # get plotmatrix for overview
    def getPlot(self):   
        plotmatrix = np.zeros([self.initComp._count, 21])
        plotmatrix[:,0] = self.initComp._tins
        plotmatrix[:,1] = self.sumAll.container[:,0]
        plotmatrix[:,2] = self.sumAll.container[:,1]
        plotmatrix[:,3] = self.sumAll.container[:,2]
        plotmatrix[:,4] = self.sumAll.container[:,3]
        plotmatrix[:,5] = self.BAT.compContain.container[:,2]
        plotmatrix[:,6] = self.BAT.compContain.container[:,0]
        plotmatrix[:,7] = self.BAT.compContain.container[:,1]
        plotmatrix[:,8] = self.SOLAR.compContain.container[:,0]
        plotmatrix[:,9] = self.SOLAR.compContain.container[:,1]
        plotmatrix[:,10] = self.SOLAR.compContain.container[:,2]
        plotmatrix[:,11] = self.PCU.compContain.container[:,0]
        plotmatrix[:,12] = self.PCU.compContain.container[:,1]
        plotmatrix[:,13] = self.OBC.compContain.container[:,0]
        plotmatrix[:,14] = self.OBC.compContain.container[:,1]
        plotmatrix[:,15] = self.AOCS.compContain.container[:,0]
        plotmatrix[:,16] = self.AOCS.compContain.container[:,1]
        plotmatrix[:,17] = self.COMUHF.compContain.container[:,0]
        plotmatrix[:,18] = self.COMUHF.compContain.container[:,1]
        plotmatrix[:,19] = self.COMSBAND.compContain.container[:,0]
        plotmatrix[:,20] = self.COMSBAND.compContain.container[:,1]
        return plotmatrix
    
    
    
        # shortent notation
        # Dissipation Power has to be performed after calcEconsum and is needed for further FEM Calc
    def getDissipation(self):
        plotmatrix = np.zeros([self.initComp._count, self.count])
        for i in range(0, self.count):
            plotmatrix[:,i] = self.dissiContain[i].container[:,1]
        return plotmatrix
    
            # Heatflux calc from Dissipation and specified Area [mm^2] but returns in [W/m^2] has to be performed after calcEconsum and is needed for further FEM Calc
    def getHeatflux(self):
        plotmatrix = np.zeros([self.initComp._count, self.count])
        for i in range(0, self.count):
            plotmatrix[:,i] = self.dissiContain[i].container[:,1]  / self.dissiContain[i].getArea() * 10.**6
        return plotmatrix
    
        # returns name for further FEM Calc (Output)
    def getName(self):
        name = []
        for i in range(0, self.count):
            name.append(self.dissiContain[i].getName())
        return name
        
    
    
              # get solar load for further heat reduction
    def getSolarload(self):
        return self.SOLAR.compContain.container[:,1]
    
            # get dissipation for further thermal simuolation
    def getSumDissipation(self):
        return self.sumAll.container[:,1]