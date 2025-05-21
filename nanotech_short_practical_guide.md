# MD for Nanotech quick 'n clean

## Clone this git repo    
To work with this code make sure to clone this exact git repo:    
`git clone <url of repo>`

## Using Python Packages
Packages include functions, objects, etc. that are not natively included in python.    
In the course of this practical we will use a number of packages i.e. `ase` for the actual simulation, `numpy` for matrix calculations, `scipy` for curve fitting. To use functions and objects we need to import them:    
`import <insert_package_name_here> as <name_of_alias>`     
e.g `import numpy as np`    
if we want to use just one particular function we can use:    
`from <package> import <function_or_object_from_package>`   
e.g. `from ase.optimize import FIRE`    

## Setting up the environment    
Before any of the scripts will run, make sure to initialize a virtual environment with all the required packages.     
To do this open the directory of this code in a terminal. If you cloned it to your home directory you should acces it via:    
`cd <name of clone dir>`     
Then set up the environment here (this one sets up an environment called env).         
`python -m venv env`    
Now activate the environment:     
`source env/bin/activate`    
After this your command line should look like thise `(env) <yourusername>@<nameofyourmachine>:/`     
To now include all the packages for the python scripts run:    
`pip install -r requirements.txt`    


## Running the code for determining the lattice constants     
Run the script:    
`python3 minimize_energy_script.py`     
The script will ask you to provide initial guesses for the lattice constants and will calculate 300 data points in a range of -0.3 Angstrom and +0.3 Angstrom of your initial guess.      
The strain energydentsity csv files are written to the folder `data_min_energy/`     
Plot the data and find the minima.
