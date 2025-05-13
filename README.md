# practical_atomistic_simulation
Practical exercises for atomistic simulation utilizing ASE
make sure to be in the folder of your project    
for setting up the environment use:    
`python3 -m venv create env`             
then to enter the environment:    
`source env/bin/activate`    
after this install packages (as listed in the requirements.txt file)    
`pip install -r requirements.txt´        

if conda is used:
to set up the environment use conda - start with  
`conda env create --name atomistic --file=environment.yml`  
`conda activate atomistic`

if work has been done adding new packages, etc. changing the environment use  
`conda env export > environment.yml`  
to update .yml

the python files are associated with different Tasks in the Problem Statement  
minimize_energy.py -> Task 3.2
elastic_constants.py -> Task 3.5
