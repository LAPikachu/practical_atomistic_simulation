# MD for Nanotech quick 'n clean

## Using Python Packages
Packages include functions, objects, etc. that are not natively included in python.    
In the course of this practical we will use a number of packages i.e. `ase` for the actual simulation, `numpy` for matrix calculations, `scipy` for curve fitting. To use functions and objects we need to import them:    
`import <insert_package_name_here> as <name_of_alias>`     
e.g `import numpy as np`    
if we want to use just one particular function we can use:    
`from <package> import <function_or_object_from_package>`   
e.g. `from ase.optimize import FIRE`    