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

directory = 'data_elastic_constants/'

#this object stores a strain tensor for the 3 different types we investigate
class StrainTensor:
    def __init__(self):
        self.tensor = np.zeros((3,3))
        
    def make_tensor(self, type, eps):
        if type == 'uniax':
            self.uniax(eps)
        elif type == 'biax':
            self.biax(eps)
        elif type == 'shear':
            self.shear(eps)
        else:
             raise ValueError(f"No deform type {type}")
    def uniax(self, eps):
        self.tensor[0,0] = eps

    def biax(self, eps):
        self.tensor[0,0] = eps
        self.tensor[1,1] = eps

    def shear(self, eps):
        self.tensor[0,1] = 2*eps
        self.tensor[1,0] = 2*eps

def deform_cell(cell, strain_tensor):
    cell_array = cell.cell[:]
    # don't use * for multiplying matrices
    cell_array_new = np.dot(cell_array,(np.identity(3)+strain_tensor))
    cell.cell[:] = cell_array_new
    return cell

if __name__ == '__main__':
    #relax initial structur
    with open("data_min_energy\lattice_consts", "r") as file: #load lattice constant
        lattice_const_dict = json.load(file) 
    lattice_const = lattice_const_dict['W_3rd_ord']['fcc']
    cell = create_cubic_lattice('fcc', lattice_const)
    cell.calc = EMT()
    #relax cell  
    relax_cell =  BFGS(cell, trajectory=f'{directory}relaxed_cell.traj')
    relax_cell.run()
    traj_relax = Trajectory(f'{directory}relaxed_cell.traj', 'r')
    cell_relaxed = traj_relax[-1]

    deform_types = ['uniax', 'biax', 'shear']
    epsilons = np.linspace(-0.02,0.02, 20) #watch out: no big vals for eps 

    #deform the cell 
    for type in deform_types:
        #this open/close dance is necessary to load and then start work on traj
        traj = Trajectory(f'{directory}{type}.traj', 'w')
        traj.write(cell_relaxed)
        traj = Trajectory(f'{directory}{type}.traj', 'r')
        cell = traj[-1]
        traj = Trajectory(f'{directory}{type}.traj', 'w')

        tensor = StrainTensor()

        for eps in epsilons:
            tensor.make_tensor(type, eps)
            cell = deform_cell(cell, tensor.tensor)
            cell.calc = EMT()
            cell.get_potential_energy()
            traj.write(cell)

    #fit and find min
    

