# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:20:07 2017
Demonstration File early alpha release
@author: Danilo Költzsch
"""
# for demonstration with open geometry libary with heat flux (works) 

import pandas as pd
import numpy as np
import os
import time

# Lotte packages used for demonstration
from orbit import orbit
from sat import satorient 
from surfaceheat2 import surfaceheatflux 
from openLibrary import *
from SimpleSolver import ThermalOneNodeSolver
from Output import OutputGenerator 




#set file for save and eventually creates folder
file_path = "D:\Python\Plot\demoB4-2\\"
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
   os.makedirs(directory)
   
# note for data identifier (not more than 14 chars)
note = 'CubeSat freetumble'
note2 =  note + ' ' + time.strftime("%H-%M %d.%m.%Y")






# init from library of geometry
OpenSat = Cube()
surfaceOpenSat  = OpenSat._sf

# init Orbit data
orb = orbit()


textISS = """
1 25544U 98067A   17046.50386891  .00016717  00000-0  10270-3 0  9005
2 25544  51.6424 284.9118 0006538 195.9975 164.0972 15.54379348  2824
"""





# orbit setup
# defines time frame for orbit
# time procedure, sets time frame for simulation in tsta in utc; tend in min; count integer    
# def timeframe(self, tend, count, year, month, day, hour, minute, sec):
# defines the duration for propagation in minutes (tend)
duration = 400.
# defines the number of timesteps for specified duration (count)
stepcount = 2001
orb.timeframe(duration, stepcount, 2016, 9, 9, 19, 0, 0)
# def Berlin ILR for contact simulation
orb.defpositionLocal('52.5145 N', '13.3237 E')
# propagation for specied TLE as text for timeframe
orb.propagation(textB4)
# caclulation if satellite is illuminated or not  
orb.illumination()


# creates satellite attitude object
sator = satorient()
# redifnes orbit by orbit object, which was previous defined
sator.defOrbit(orb)


# defines satellite local flight vectors for nadir mode (can be modified)
nadir = np.array([0,0,1])
fly =   np.array([1,0,0])

# actual possible modes are:
#  - simplenadir (nadir pointing flight fix)
#  - nadirspin (nadir pointing constant spin)
#  - freetumble (constant spin)
#  - sunspin (sun pointing with spin)
#sator.defNadirV(nadir, fly)

# calculates simple nadir mode
#sator.simplenadir()


# example for free spin:
# second alternative mode (overwrites nadir mode, just for deomstration)    
# tumble rate in °/s about x y z  (alpha beta gamma)
tumblerateX = 0.004
tumblerateY = 0.001
tumblerateZ = 0.02

sator.defFreetumbleV(nadir, fly,tumblerateX,tumblerateY,tumblerateZ)
sator.freetumble()






# inits heatflux object for environment simulation based on stem
satheat = surfaceheatflux()
# specifies environmentdata for solar , earth infrared and albedo (here medium, possible cold/hot)
satheat.setEnvironmentByCase('medium', 'medium',  'medium')
# calculates outer heat fluxes in w/m^2 by orbit, sat attitude and surface 
satheat.surfaceHeatflux(orb, sator, surfaceOpenSat)

# advanced Dissipation with Solar Charging based on illumination simulation performed before
dissi = StandardChargeElectroSystem(satheat)

# the definition is modus based (defined in Library) and is controlled here
count = dissi.initComp._count
c1 = count / 2
for i in range(0, count):
    if i < c1 and orb.getDistance(i) < 1500.:
        dissi.calcEconsum(i, 1)
    elif i < c1:
        dissi.calcEconsum(i, 0)
    elif i >= c1 and orb.getDistance(i) < 1500.: 
        dissi.calcEconsum(i, 3)     
    else:
        dissi.calcEconsum(i, 2)
        
dissi_overview = dissi.getPlot()

# reduces solar heatfluxes on solar arrays depending on the sat systems powercosumption and the load on the solar cells
satheat.reduceHeatflux(dissi.getSolarload())


# simple thermal node solver for middle temperature (needed for pre processing and object comparison)
tempcalc = ThermalOneNodeSolver()
# defines the satellite object (geometry and mass)
tempcalc.defParam(surfaceOpenSat, 1.1)
# calculation with reduced heatfluxes from objects satheat and dissipation sum and starting temperature
tempcalc.calcTempRed(satheat, dissi.getSumDissipation(), 0)
# reads the temperature solution
tempsolu= tempcalc._container


# defines simplified output from objects satheat and dissipation 
ouroutput = OutputGenerator(satheat, dissi, note)
# defines dissipation case (output style)
# returnes Dissipation defined as Power in [W] [0] or as heatflux related to the surface size [W/m^2]
ouroutput.dissipationCase(0)
# defines heatflux case (output style)
# writes heatflux with case [ [0] standard heatflux [W/m^2], [1] reduced heatflux [W/m^2],  [2] standard power (related to defined surface [W],  [3] reduced power (related to defined surface [W]; [4] standard heatflux albedo respecting [W/m^2], [5] reduced heatflux albedo respecting [W/m^2]]
ouroutput.heatfluxCase(1)
# possible data reduction defined by number of intended timesteps
# reduction for all relevant components for less timesteps    
ouroutput.reduction(200)
# defines typ of outputfile here excel sheet stored at filepath with name 'demo_output'
# possible other forms
# - writeCSV
# - writeANSYS
ouroutput.writeXLSX('demo_output',file_path)


# Export into XLS of Orbit data

light = orb.getIllumination()
cols = ['time [s]', 'sun sat vector x [km]', 'sun sat vector y [km]', 'sun sat vector z [km]', 'sat vector x [km]', 'sat vector y [km]', 'sat vector z [km]','sun vector x [km]', 'sun vector y [km]', 'sun vector z [km]', 'vsat vector x [km/s]', 'vsat vector y [km/s]', 'vsat vector z [km/s]', 'distance [km]', 'illumination [0/1]', ]
df = pd.DataFrame(data=light, columns=cols)
writer = pd.ExcelWriter(file_path+'orbit_illu.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name=note2)
writer.save()

# Export into XLS of Satellite Attitude data

satdcm = sator.getSatDCM()
cols2 = ['time [s]', 'xx vector', 'xy vector', 'xz vector', 'yx vector', 'yy vector', 'yz vector','zx vector', 'zy vector', 'zz vector']
df2 = pd.DataFrame(data=satdcm, columns=cols2)
writer2 = pd.ExcelWriter(file_path+'sat_dcm.xlsx', engine='xlsxwriter')
df2.to_excel(writer2, sheet_name=note2)
writer2.save()

# Export into XLS of Satellite Middle Temperature

cols3 = ['Time in s','Temperature in C','Nettowork in W','Dissipation in W','Incomming Radiation in W','Outgoing Radiation in W']
df3 = pd.DataFrame(data=tempsolu, columns=cols3)
writer3 = pd.ExcelWriter(file_path+'temp_calc_red.xlsx', engine='xlsxwriter')
df3.to_excel(writer3, sheet_name=note2)
writer3.save()


# Export into XLS of Satellite Esystem
cols4 = ['Time in s','Econsumption All in W','Dissipation All in W','Econsumption All without Charge in W','Mode','Charge in Bat in Ws','Econsumption Bat in W','Dissipation Bat in W','Egeneration Solar in W','SolarEnergy Consumption in W','Incomming Solarpower [incl Alb.] in W','Econsumption PCU in W','Dissipation PCU in W','Econsumption OBC in W','Dissipation OBC in W','Econsumption AOCS in W','Dissipation AOCS in W','Econsumption COMUHF in W','Dissipation COMUHF in W','Econsumption COMSBAND in W','Dissipation COMSBAND in W']
df4 = pd.DataFrame(data=dissi_overview, columns=cols4)
writer4 = pd.ExcelWriter(file_path+'OpenSat_esystem.xlsx', engine='xlsxwriter')
df4.to_excel(writer4, sheet_name=note2)
writer4.save()



