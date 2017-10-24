# LOTTE
"Low earth Orbit Tool for Thermal Environment" - Thermal Space Simulation Tool Chain for Finite Element Analysis boundary conditions in Python

soon available

1. please add the following packages to Python:
 - pip install numpy
 - pip install scipy
 - pip install pandas
 - pip install skyfield

- for windows shapely:
  - conda config --add channels conda-forge
  - conda install shapely
  - or direct wheels installation  http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely  see also https://pypi.python.org/pypi/Shapely

 - for other os shapely
   - pip install shapely

2. Discription will be soon available

 - Use demonstration for introduction (well filled with comments)

3. Already ready:

- [x] openLibrary - simple geometry and dissipation systems for small cubesats
- [x] quateul - Quaternion calc in Python
- [x] orbit - Orbit propagation with skyfield
- [x] sat   - Sat Attitude Propagation
- [x] heatflux2 - works for simple geometries with solar reduction and geometric shader implemented yet
- [x] GeometricShader - simple geo shaeders partly using shapely
- [x] dissipation - works with battery charging and solar energy generation, standard charge system is the advanced electrical simulation yet
- [x] SimpleSolver - solution to middle temperatures from generated solution by one node model for better initial temperature and extrem - [x] orbit detection
- [x] output - export data in different writing formats (heatflux and power) also reduces exported data to an intendend amount of steps, exports into csv, xlsx and ASNYS readable file 

- [x] openlibrary - including surface description and dissipation for small sats (4 satellites available 1u 2u deployed octa)
