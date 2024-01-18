import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read
from minimize_energy import create_cubic_lattice
from ase.optimize import BFGS

#print(cell.cell[:]) # cell.cell only gives a 1-dim array type output :/ ?

if __name__ == '__main__':
    #make cell
    with open("data_min_energy\lattice_consts", "r") as file: #load lattice constant
        lattice_const_dict = json.load(file) 
    lattice_const = lattice_const_dict['W_3rd_ord']['fcc']
    cell = create_cubic_lattice('fcc', lattice_const)
    cell.calc = EMT()
    #relax cell  
    relax_cell =  BFGS(cell, trajectory='data_elastic_constants/relaxed_cell.traj')
    relax_cell.run()

