# LOTTE
"Low earth Orbit Tool for Thermal Environment" - Thermal Space Simulation Tool Chain for FEA boundary conditions in Python

please add the following packages to Python:
pip install numpy
pip install scipy
pip install pandas
pip install skyfield

for windows:
conda config --add channels conda-forge
conda install shapely

or direct wheels installation  http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely  see also https://pypi.python.org/pypi/Shapely

for other
pip install shapely

Discription will be soon available

Use demonstration for introduction (well filled with comments)

Already ready:

openLibrary - simple geometry and dissipation systems for small cubesats
quateul - Quaternion calc in Python
orbit - Orbit propagation with skyfield
sat   - Sat Attitude Propagation
heatflux2 - works for simple geometries with solar reduction and geometric shader implemented yet
GeometricShader - simple geo shaeders partly using shapely
dissipation - works with battery charging and solar energy generation, standard charge system is the advanced electrical simulation yet
SimpleSolver - solution to middle temperatures from generated solution by one node model for better initial temperature and extrem orbit detection
output - export data in different writing formats (heatflux and power) also reduces exported data to an intendend amount of steps, exports into csv, xlsx and ASNYS readable file 

openlibrary - including surface description and dissipation for small sats (4 satellites available 1u 2u deployed octa)
