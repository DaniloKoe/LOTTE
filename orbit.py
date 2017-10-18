# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 14:01:44 2017
Illumination Test for Orbit Simulation
@author: Danilo Költzsch
"""
import math 
import numpy as np
from skyfield.positionlib import ICRF
from skyfield.api import Star, load, Topos, EarthSatellite
import decimal as dec


class orbit:

  
    
    # init procedur
   def __init__(self):
        self._ts = load.timescale()
        self._count = 2
        self._year = 2000
        self._month = 1
        self._tincr = 10
        self._utcday = 1
        self._tutcday = np.array([self._utcday,self._utcday + self._tincr / 86400])
        self._tutcday2 = np.array([self._utcday + 0.1 / 86400,self._utcday + (self._tincr + 0.1) / 86400])
        self._tutc = self._ts.utc(self._year, self._month, self._tutcday)
        self._tutc2 =self._ts.utc(self._year, self._month, self._tutcday2)
        self._tins = np.zeros((self._count))
        # Orbit Stuff
        self._TLE = """
ISS (ZARYA)
1 25544U 98067A   17046.50386891  .00016717  00000-0  10270-3 0  9005
2 25544  51.6424 284.9118 0006538 195.9975 164.0972 15.54379348  2824
"""  
        # gps earth radius in km
        self._rearth =  	6378.137
        # initial skyfield prop
        self._planets = load('de421.bsp')
        self._earth = self._planets['earth']
        self._earthpos  = self._earth.at(self._tutc)
        self._earthpos2 = self._earth.at(self._tutc2)
        self._line1, self._line2 = self._TLE.splitlines()[-2:]
        self._sat = self._earth + EarthSatellite(self._line1, self._line2, name="TUB Sat")
        self._sunv = np.transpose(self._earthpos.position.km)
        self._sunv2 = np.transpose(self._earthpos2.position.km)
        self._satobj =  self._sat.at(self._tutc)
        self._satobj2 = self._sat.at(self._tutc2)
        self._sunsatv  = np.transpose(self._satobj.position.km)
        self._sunsatv2 = np.transpose(self._satobj2.position.km)
        self._satv  = np.subtract(self._sunv,self._sunsatv)
        self._satv2 = np.subtract(self._sunv2,self._sunsatv2)
        self._satvel = np.subtract(self._satv2,self._satv[:,:])*10      
        self._light = np.zeros((self._count))
        self._mm  = 15.54379348 / 24 / 60 / 60
        self._T   =  1 / self._mm
        self._muEarth = 3.986004418 * 10**14
        # Berlin as EarthPosition
        self._earthLoc = self._earth + Topos('52.5145 N', '13.3237 E')
        self._earthLocPos = self._earthLoc.at(self._tutc)
        self._earthLocPosv= np.transpose(self._earthLocPos.position.km)
        self._earthLocPosSatv = np.subtract(self._earthLocPosv, self._sunsatv)
        self._earthLocDistance = (self._earthLocPos - self._satobj).distance().km
        #print('Orbit Object initialised')
    
       

        
   # time procedure, sets time frame for simulation in tsta in utc; tend in min; count integer    
   def timeframe(self, tend, count, year, month, day, hour, minute, sec):
        self._count = count
        # declare arrays new before flooded with data
        self._tins = np.zeros(self._count)
        self._light = np.zeros(self._count)
        self._tutcday = np.zeros(self._count)
        self._tutcday2 = np.zeros(self._count)
        # utc format in days
        self._utcday = day + hour / 24 + minute / 1440 + sec / 86400 
        # set new time from specification
        self._tincr = tend / (self._count - 1) * 60
        self._month = month
        self._year =  year
        for x in range(0, self._count):
            self._tutcday[x] = self._utcday + x * self._tincr / 86400
            # smallvector for speedvector
            self._tutcday2[x] = self._utcday  + (x  * self._tincr + 0.1) / 86400
            self._tins[x] = x * self._tincr
        self._tutc = self._ts.utc(self._year, self._month, self._tutcday)
        self._tutc2 = self._ts.utc(self._year, self._month, self._tutcday2)
        print('Time modified as specified')

    
   def propagation(self, text):
         # produces Satellite with TLE's zext
        self._earthpos = self._earth.at(self._tutc)  
        self._earthpos2 = self._earth.at(self._tutc2) 
        self._line1, self._line2 = text.splitlines()[-2:]
        self._sat = self._earth + EarthSatellite(self._line1, self._line2)
        # From the center of the Solar System (Barycentric)
        tempsunvt  = self._earthpos.position.km
        tempsunvt2  = self._earthpos2.position.km
        self._sunv = np.transpose(tempsunvt) * -1.
        self._sunv2 = np.transpose(tempsunvt2) * -1.
        # From the center of the Earth (Geocentric)
        self._satobj  = self._sat.at(self._tutc)
        self._satobj2 = self._sat.at(self._tutc2)
        tempsunsatvt     = self._satobj.position.km
        tempsunsatv2t    = self._satobj2.position.km
        self._sunsatv    = np.transpose(tempsunsatvt) * -1.
        self._sunsatv2   = np.transpose(tempsunsatv2t) * -1.
        self._satv  = np.subtract(self._sunv,self._sunsatv)
        self._satv2 = np.subtract(self._sunv2,self._sunsatv2)
        self._earthLocPos = self._earthLoc.at(self._tutc)
        self._earthLocPosv= np.transpose(self._earthLocPos.position.km)
        self._earthLocPosSatv = np.subtract(self._earthLocPosv, self._sunsatv)
        self._earthLocDistance = (self._earthLocPos - self._satobj).distance().km
        # Generates velocity vector
        self._satvel = np.subtract(self._satv2,self._satv)*10 
        print('Orbit propagated')
        

        # propagation must be performed before
   def getPropagation(self):             
        plotmatrix = np.zeros([self._count,11])
        plotmatrix[:,0]=self._tins[:]
        plotmatrix[:,1]=self._sunsatv[:,0]
        plotmatrix[:,2]=self._sunsatv[:,1]
        plotmatrix[:,3]=self._sunsatv[:,2]
        plotmatrix[:,4]=self._satv[:,0]
        plotmatrix[:,5]=self._satv[:,1]
        plotmatrix[:,6]=self._satv[:,2]
        plotmatrix[:,7]=self._satvel[:,0]
        plotmatrix[:,8]=self._satvel[:,1]
        plotmatrix[:,9]=self._satvel[:,2]
        plotmatrix[:,10]=self._earthLocDistance[:]
        return plotmatrix
        
    
    
        
    # Illumination test procedure, checks with sun and satellite vector in km for illumination
   def illumination(self):
        # satshape has to be the same as sunshape
        # test for debugging
        #print(light)
        for i in range(0, self._count):
            #  generates a therm of cubic equation
            a = self._sunsatv[i,0]**2 + self._sunsatv[i,1]**2  + self._sunsatv[i,2]**2
            #  generates b therm of cubic equation
            b = - 2 * ( self._sunsatv[i,0] * self._sunv[i,0] +  self._sunsatv[i,1] * self._sunv[i,1] + self._sunsatv[i,2] * self._sunv[i,2])
            #  generates c therm of cubic equation
            c = self._sunv[i,0]**2 + self._sunv[i,1]**2 + self._sunv[i,2]**2 - self._rearth**2
            #  generates d therm for simplification  
            d = (b * b - 4 * a * c )  
            if d > 0:
                s = (- b - math.sqrt(d)) /(2 * a)
                if s > 1:
                    # satellites is located in rear of earth, not illuminated
                    self._light[i] = 0
                else:    
                    # satellites is located in front of earth, illuminated
                    self._light[i] = 1
            elif d == 0:
                # satellite can not be on earth surface
                self._light[i] = 2
            elif d < 0:
                # satellites is located beside earth, illuminated
                self._light[i] = 1
            else:
                # no possible solution
                self._light[i] = 3
        print('Illumination propagated')    
        

        
        
         # illuminationn must be performed before
   def getIllumination(self):     
        plotmatrix = np.zeros([self._count,15])
        plotmatrix[:,0]=self._tins[:]
        plotmatrix[:,1]=self._sunsatv[:,0]
        plotmatrix[:,2]=self._sunsatv[:,1]
        plotmatrix[:,3]=self._sunsatv[:,2]
        plotmatrix[:,4]=self._satv[:,0]
        plotmatrix[:,5]=self._satv[:,1]
        plotmatrix[:,6]=self._satv[:,2]
        plotmatrix[:,7]=self._sunv[:,0]
        plotmatrix[:,8]=self._sunv[:,1]
        plotmatrix[:,9]=self._sunv[:,2]
        plotmatrix[:,10]=self._satvel[:,0]
        plotmatrix[:,11]=self._satvel[:,1]
        plotmatrix[:,12]=self._satvel[:,2]
        plotmatrix[:,13]=self._earthLocDistance[:]
        plotmatrix[:,14]=self._light[:]
        return plotmatrix
   
    # TLE Generator by Thomas Meschede
   def buildtle(self, inc, lan, ecc, aop, ma, mm, yr, day):
        # revolutions per sec
        self._mm  = mm / 24 / 60 / 60
        # revolution Time 
        self._T   = 1 / self._mm
        """inc: inclination (degrees)
           lan: longitude of ascending node (degrees, also right ascension of the ascending node)
           ecc: eccentricity (decimal point assumed)
           aop: argument of perigee (degrees)
           ma:  mean anomaly (degrees)
           mm:  mean motion (revolutions per day)"""
        tmpl="""
1 25544U 98067A   {yr:d}{day:012.8f} -.00002182  00000-0 -11606-4 0 2927
2 99999  {inc:07.4f} {lan:08.4f} {ecc:07d} {aop:08.4f} {ma:08.4f} {mm:11.8f}563537
""".format(inc=inc, lan=lan, ecc=int(ecc*10e6), aop=aop, ma=ma, mm=mm, yr=yr, day=day)
        self._TLE = tmpl
        return tmpl  
    
    
    #  orbit time per semi-major axis
   def calcT(self, sma):
        self._T = 2 * np.pi * np.linalg.sqrt(sma**3 / self.__muEarth)
        self._mm = 1 / self._T
        return self._T
   
    # extracts Mean Motion from interior TLE  (mm in 1 / s)
   def extractMeanMotion(self):
        tlelocV = self._TLE.split()
        meanmotionstr = tlelocV[-2]
        self._mm = dec.Decimal(meanmotionstr) / 24 / 60 / 60
        self._T = 1 / self._mm
        print('Mean Motion from intern TLE extracted: (1/s)', self._mm)
        
       # calcs radius (flight hight)
   def getRadius(self, step):
        h =  np.linalg.norm(self._satv[step,:]) - self._rearth
        return h
   
    # extracts Mean Motion from given Text TLE and returns it (mm in 1 / day)
   def dirExtractMeanMotion(self, text):
        tlelocV = text.split()
        meanmotionstr = tlelocV[-2]
        meanmotion = dec.Decimal(meanmotionstr) / 24 / 60 / 60
        print('Mean Motion from Text TLE extracted: (1/day)', meanmotion)
        return meanmotion
    
    # defines earth local position with length and width
   def defpositionLocal(self, textLength, textWidth):
        self._earthLoc = self._earth + Topos(textLength, textWidth)
        self._earthPos = self._earthLoc.at(self._tutc)
        self._earthPosv= np.transpose(self._earthPos.position.km)
    
    # returns distance between earth local position and satellite for time step
   def getDistance(self, step):
       return self._earthLocDistance[step]
        
        

        
#if __name__ == "__main__":
#    candle =  orbit()
        