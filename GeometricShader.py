# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:55:05 2017
Simple Geometric Shader
@author: Danilo Költzsch
"""

# for all options needs python 3.5

import numpy as np
#from geosat import SurfaceGeometry as sfg

# new shading possibilities with shapely, shapely runnings currently only under python 3.5
from shapely.geometry import Polygon

    # dummyShader
class dummyShader:
    
       # inits solver by definig geometry
    def __init__(self):
        self.length         = 0
        self.height         = 0
        #print('Dummy Shader inited')
      
     # in -1 direction 1 
    def calcShader(self, angle1, angle2):
        factor = 1.
        return factor



    # one wall on one surface
class oneSideShader:
    
       # inits solver by definig geometry # in -1 direction 1
    def __init__(self, length1, height1):
        self.length         = length1
        self.height         = height1
        #print('One Side Shader inited')



       
     # in -1 direction 1 
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2): 
            shadedlength = self.height * np.tan(-angle1)
            if shadedlength <= self.length and shadedlength > 0:
                factor = (self.length - shadedlength) / self.length
            elif shadedlength > self.length:
                factor = 0
            else:
                factor = 1.
        elif angle1 < (np.pi / 2):
            factor = 1. 
        else:
            factor = 0
        return factor
    
    # two walls beside on surface
class twoSideShader:

           # inits solver by definig geometry # in -1 direction 1 and -1 direction 2
    def __init__(self, length1, length2, height1, height2):
        self.length1         = length1
        self.height1         = height1
        self.length2         = length2
        self.height2         = height2
        #print('Two Side Shader inited')
       
        
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2): 
            shadedlength = self.height1 * np.tan(-angle1)
            if shadedlength <= self.length1 and shadedlength > 0:
                factor1 = (self.length1 - shadedlength) / self.length1
            elif shadedlength > self.length1:
                factor1 = 0
            else:
                factor1 = 1.
        elif angle1 < (np.pi / 2):
            factor1 = 1.
        else:
            factor1 = 0
        if angle2 < 0 and angle2 > (- np.pi / 2): 
            shadedlength = self.height2 * np.tan(-angle2)
            if shadedlength <= self.length2 and shadedlength > 0:
                factor2 = (self.length2 - shadedlength) / self.length2
            elif shadedlength > self.length2:
                factor2 = 0
            else:
                factor2 = 1.
        elif angle2 < (np.pi / 2):
            factor2 = 1.
        else:
            factor2 = 0
        factor = factor1 * factor2        
        return factor
    
    # enclosed surfaces by 4 walls
class enclosedShader:

    
       # inits solver by definig geometry # in -1 direction 1 and -1 direction 2
    def __init__(self, length1, length2, height1, height2, height3, height4):
        self.length1         = length1
        self.height1         = height1
        self.length2         = length2
        self.height2         = height2
        self.height3         = height3
        self.height4         = height4        
        #print('Enclosed Shader inited')


       
        
    def calcShader(self, angle1, angle2):
        if angle1 <= 0 and angle1 > (- np.pi / 2): 
            shadedlength = self.height1 * np.tan(-angle1)
            if shadedlength <= self.length1 and shadedlength > 0:
                factor1 = (self.length1 - shadedlength) / self.length1
            elif shadedlength > self.length1:
                factor1 = 0
            else:
                factor1 = 1.
        elif angle1 > 0 and angle1 < (np.pi / 2):
            shadedlength = self.height3 * np.tan(angle1)
            if shadedlength <= self.length1 and shadedlength > 0:
                factor1 = (self.length1 - shadedlength) / self.length1
            elif shadedlength > self.length1:
                factor1 = 0
            else:
                factor1 = 1. 
        else:
            factor1 = 0
        if angle2 <= 0 and angle2 > (- np.pi / 2): 
            shadedlength = self.height2 * np.tan(-angle2)
            if shadedlength <= self.length2 and shadedlength > 0:
                factor2 = (self.length2 - shadedlength) / self.length2
            elif shadedlength > self.length2:
                factor2 = 0
            else:
                factor2 = 1.
        elif angle2 > 0 and angle2 < (np.pi / 2):
            shadedlength = self.height4 * np.tan(angle2)
            if shadedlength <= self.length2 and shadedlength > 0:
                factor2 = (self.length2 - shadedlength) / self.length2
            elif shadedlength > self.length2:
                factor2 = 0
            else:
                factor2 = 1. 
        else:
            factor2 = 0
        factor = factor1 * factor2        
        return factor
    
    
    # simple shader for small wall fragments relatad to shaded surface
class smallObjectShader:
        
       # inits solver by definig geometry
    def __init__(self, length1, height1, areawidth, objectwidth):
        self.length         = length1
        self.height         = height1
        self.objectwidth    = objectwidth
        self.areawidth      = areawidth        
        #print('Small Object Shader inited')


       
     # in -1 direction 1 
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2): 
            shadedlength = self.height * np.tan(-angle1)
            if shadedlength <= self.length and shadedlength > 0:
                factor = 1 - (self.objectwidth * shadedlength) / (self.areawidth * self.length)
            elif shadedlength > self.length:
                factor = 1 - self.objectwidth / self.areawidth
            else:
                factor = 1.
        elif angle1 < (np.pi / 2):
            factor = 1. 
        else:
            factor = 0
        return factor
    
    
    # object shader for small wall fragments relatad to shaded surface
class ObjectShader:
        
       # inits solver by definig geometry
    def __init__(self, length, height, areawidth, objectwidth, posobject):
        self.length         = length
        self.height         = height
        self.objectwidth    = objectwidth
        self.areawidth      = areawidth
        self.posobject      = posobject        
        #print('Object Shader inited')

    


       
     # using shapely for cool geometrical shader v3.5 needed currently
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2) and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            # geometrical calc of polygon edges
            shadedlength = self.height * np.tan(-angle1)
            shadededge1   = self.posobject + shadedlength * np.tan(angle2)
            shadededge2   = shadededge1 + self.objectwidth
            posobject2    = self.posobject + self.objectwidth
            # shapley magic
            shadow = Polygon([(self.posobject,0),(posobject2,0),(shadededge2,shadedlength),(shadededge1,shadedlength)])
            surface = Polygon([(0,0),(self.areawidth,0),(self.areawidth,self.length),(0,self.length)])
            shadedsurface = surface.intersection(shadow)
            if surface.intersects(shadow):
                factor = 1. - shadedsurface.area / surface.area 
                if factor < 0:
                    factor = 0.
            else:
                factor = 1.
        elif angle1 > 0 and angle1 < (np.pi / 2)  and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            factor = 1. 
        else:
            factor = 0
        return factor



    # wall shader, facing a object and satellite one side (-) / shader  relatad to shaded surface
class facingObjectShader:
        
       # inits solver by definig geometry
    def __init__(self, length, height, areawidth, objectwidth, heightobject, posobject, distanceobject):
        self.length         = length
        self.height         = height
        self.objectwidth    = objectwidth
        self.areawidth      = areawidth
        self.heightobject   = heightobject
        self.posobject      = posobject
        self.disobject      = distanceobject        
        #print('Inner Wall Shader inited')

       
     # using shapely for cool geometrical shader v3.5 needed currently
    def calcShader(self, angle1, angle2):
        surface = Polygon([(0,0),(self.areawidth,0),(self.areawidth,self.length),(0,self.length)])
        if angle1 < 0 and angle1 > (- np.pi / 2): 
            # geometrical calc of polygon edges
            shadedlength = self.height * np.tan(-angle1)
            if shadedlength <= self.length and shadedlength > 0:
                shadowwall = Polygon([(0,0),(self.areawidth,0),(self.areawidth,shadedlength),(0,shadedlength)])
                objectpos1 = self.disobject  * np.tan(-angle1)
                objectpos2 = self.disobject  * np.tan(angle2)
                objp21     = objectpos2 + self.posobject
                objp22     = objp21 + self.objectwidth
                objp11     = objectpos1
                objp12     = objectpos1 + self.heightobject
                shadowobj  = Polygon([(objp21,objp11),(objp22,objp11),(objp22,objp12),(objp21,objp12)])
                shadow     = shadowwall.union(shadowobj)
                shadedsurface = surface.intersection(shadow)
                factor = 1 - shadedsurface.area / surface.area
                if factor < 0:
                    factor = 0.
            elif shadedlength > self.length:
                factor = 0
        elif angle1 < (np.pi / 2):
            objectpos1 = self.disobject  * np.tan(-angle1)
            objectpos2 = self.disobject  * np.tan(angle2)
            objp21     = objectpos2 + self.posobject
            objp22     = objp21 + self.objectwidth
            objp11     = objectpos1
            objp12     = objectpos1 + self.heightobject
            shadow     = Polygon([(objp21,objp11),(objp22,objp11),(objp22,objp12),(objp21,objp12)])
            shadedsurface = surface.intersection(shadow)
            factor = 1 -  shadedsurface.area / surface.area
            if factor < 0:
                    factor = 0.
        else:
            factor = 0
        return factor



    
    # wall shader, facing a wall on every site and satellite one side shader  relatad to shaded surface
class innerWallShader:
        
       # inits solver by definig geometry
    def __init__(self, length, height, areawidth, facingdistance):
        self.length         = length
        self.height         = height
        self.areawidth      = areawidth
        self.distance       = facingdistance        
        #print('Inner Wall Shader inited')


       
     # using shapely for cool geometrical shader v3.5 needed currently
    def calcShader(self, angle1, angle2):
        surface = Polygon([(0,0),(self.areawidth,0),(self.areawidth,self.length),(0,self.length)])
        if angle1 > 0 and angle1 < (np.pi / 2): 
            if angle2 > 0 and angle2 < (np.pi / 2):
                # geometrical calc of polygon edges
                shaded11 = self.distance * np.tan(-angle1) 
                shaded21 = shaded11 + self.height 
                shaded12 = self.distance * np.tan(angle2)
                shaded22 = shaded12 + self.areawidth
                shadow1 = Polygon([(0,0),(0,self.height),(shaded12,shaded21),(shaded12,shaded11)])
                shadow2 = Polygon([(shaded12,shaded21),(shaded22,shaded21),(shaded22,shaded11),(shaded12,shaded11)])
                shadow     = shadow1.union(shadow2)
                shadedsurface = surface.intersection(shadow)
                factor = 1 - shadedsurface.area / surface.area
                #hier weiter
            elif angle2 == 0:
                # geometrical calc of polygon edges
                shaded11 = self.distance * np.tan(-angle1) 
                shaded21 = shaded11 + self.height 
                shadow   = Polygon([(0,shaded21),(self.areawidth,shaded21),(self.areawidth,shaded11),(0,shaded11)])
                shadedsurface = surface.intersection(shadow)
                factor = 1 - shadedsurface.area / surface.area
               #hier weiter
            elif angle2 < 0 and angle2 > (-np.pi / 2):
                # geometrical calc of polygon edges
                shaded11 = self.distance * np.tan(-angle1) 
                shaded21 = shaded11 + self.height 
                shaded12 = self.distance * np.tan(angle2)
                shaded22 = shaded12 + self.areawidth
                shadow1 = Polygon([(self.areawidth,0),(self.areawidth,self.height),(shaded22,shaded21),(shaded22,shaded11)])
                shadow2 = Polygon([(shaded12,shaded21),(shaded22,shaded21),(shaded22,shaded11),(shaded12,shaded11)])
                shadow     = shadow1.union(shadow2)
                shadedsurface = surface.intersection(shadow)
                factor = 1 - shadedsurface.area / surface.area
            else:
                factor = 0
            if factor < 0:
                factor = 0
                
        else:
            factor = 0
        return factor
    
    
        # object shader for small trapez wall fragments relatad to shaded surface
class ObjectShaderTrapez:
        
       # inits solver by definig geometry
    def __init__(self, length, height, areawidth1, areawidth2, objectwidth, posobject):
        self.length         = length
        self.height         = height
        self.objectwidth    = objectwidth
        self.areawidth1     = areawidth1
        self.areawidth2     = areawidth2
        self.areawidth21    = - ( self.areawidth2 - self.areawidth1 ) / 2
        self.areawidth22    = self.areawidth2 + self.areawidth21
        self.posobject      = posobject        
        #print('Object Shader inited')

    


       
     # using shapely for cool geometrical shader v3.5 needed currently
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2) and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            # geometrical calc of polygon edges
            shadedlength = self.height * np.tan(-angle1)
            shadededge1   = self.posobject + shadedlength * np.tan(angle2)
            shadededge2   = shadededge1 + self.objectwidth
            posobject2    = self.posobject + self.objectwidth
            # shapley magic
            shadow = Polygon([(self.posobject,0),(posobject2,0),(shadededge2,shadedlength),(shadededge1,shadedlength)])
            surface = Polygon([(0,0),(self.areawidth1,0),(self.areawidth22,self.length),(self.areawidth21,self.length)])
            shadedsurface = surface.intersection(shadow)
            if surface.intersects(shadow):
                factor = 1. - shadedsurface.area / surface.area 
                if factor < 0:
                    factor = 0.
            else:
                factor = 1.
        elif angle1 > 0 and angle1 < (np.pi / 2)  and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            factor = 1. 
        else:
            factor = 0
        return factor
    
            # object trapez shader for small wall fragments relatad to shaded surface
class ObjectTrapezShader:
        
       # inits solver by definig geometry
    def __init__(self, length, height, areawidth, objectwidth1, objectwidth2, posobject):
        self.length         = length
        self.height         = height
        self.objectwidth1   = objectwidth1
        self.objectwidth2   = objectwidth2
        self.areawidth      = areawidth
        self.objectwidth21  = - ( self.objectwidth2 - self.objectwidth1 ) / 2
        self.objectwidth22  = self.objectwidth2 + self.objectwidth21
        self.posobject      = posobject        
        #print('Object Shader inited')

    


       
     # using shapely for cool geometrical shader v3.5 needed currently
    def calcShader(self, angle1, angle2):
        if angle1 < 0 and angle1 > (- np.pi / 2) and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            # geometrical calc of polygon edges
            shadedlength = self.height * np.tan(-angle1)
            shadededge1   = self.posobject + self.objectwidth21 + shadedlength * np.tan(angle2)
            shadededge2   = shadededge1 + self.objectwidth22
            posobject2    = self.posobject + self.objectwidth1
            # shapley magic
            shadow = Polygon([(self.posobject,0),(posobject2,0),(shadededge2,shadedlength),(shadededge1,shadedlength)])
            surface = Polygon([(0,0),(self.areawidth,0),(self.areawidth,self.length),(0,self.length)])
            shadedsurface = surface.intersection(shadow)
            if surface.intersects(shadow):
                factor = 1. - shadedsurface.area / surface.area 
                if factor < 0:
                    factor = 0.
            else:
                factor = 1.
        elif angle1 > 0 and angle1 < (np.pi / 2)  and angle2 > (-np.pi / 2)  and angle2 < (np.pi / 2):
            factor = 1. 
        else:
            factor = 0
        return factor

    

       

    