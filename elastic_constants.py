import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read
from minimize_energy import create_cubic_lattice
from ase.optimize import FIRE

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
        self.tensor[0,1] = eps
        self.tensor[1,0] = eps

def deform_cell(cell_array, strain_tensor):
    # don't use * for multiplying matrices
    cell_array_new = np.dot(cell_array,(np.identity(3)+strain_tensor))
    cell.set_cell(cell_array_new, scale_atoms=True)
    return cell

if __name__ == '__main__':
    #relax initial structur
    with open("data_min_energy\lattice_consts", "r") as file: #load lattice constant
        lattice_const_dict = json.load(file) 
    lattice_const = lattice_const_dict['W_3rd_ord']['fcc']
   
    cell = create_cubic_lattice('fcc', lattice_const)
    cell.calc = EMT()
    #relax cell  
    relax_cell =  FIRE(cell, trajectory=f'{directory}relaxed_cell.traj')
    relax_cell.run()
    traj_relax = Trajectory(f'{directory}relaxed_cell.traj', 'r')
    cell_relaxed = traj_relax[-1]

    traj_check = Trajectory('check.traj', 'w')
    traj_check.write(cell_relaxed)
    
    #get initial state of relaxed cell to work with
    cell_array = cell_relaxed.cell[:]
    atoms_positions = cell_relaxed.positions
    volume = cell_relaxed.get_volume()

    deform_types = ['uniax', 'biax', 'shear']
    epsilons = np.linspace(-0.06,0.06, 20) #watch out: no big vals for eps 
    '''
    #deform the cell 
    for type in deform_types:
        cell_relaxed = read(f'{directory}relaxed_cell.traj@0', 'r')
        #traj = Trajectory(f'{directory}{type}.traj', 'r')
        #cell = traj[0]
        tensor = StrainTensor()
        traj = Trajectory(f'{directory}{type}.traj', 'w')
        
        for eps in epsilons:
            tensor.make_tensor(type, eps)
            cell = deform_cell(cell_array, tensor.tensor)
            relax_cell.run()
            cell.calc = EMT()
            cell.get_potential_energy()
            traj.write(cell)
    '''
    def fit_func(eps, p_0, p_1, p_2, p_3):
        return p_0 + p_1*eps + p_2*eps**2 + p_3*eps**3
    
    def derivative_fit_func(eps, p_1, p_2, p_3):
        return p_1 + 2*p_2*eps + 3*p_3*eps**2

    def second_derivative_fit_func(eps,p_2,p_3):
        return 2*p_2 + 6*p_3*eps

    #fit and find min
    fig1, ax1 = plt.subplots()
    fig2, (ax2, ax3) = plt.subplots(1,2)
    ax1.set_xlabel('Strain $\epsilon$')
    ax1.set_ylabel('Potential energy eV/Angstrom')
    ax2.set_ylabel('$dEpot/d\epsilon$ eV/Angstrom')
    ax2.set_ylabel('$d^{2}Epot/d\epsilon^{2}$ V/Angstrom')
    ax2.set_xlabel('Strain $\epsilon$')
    ax3.set_xlabel('Strain $\epsilon$')
    elastic_const_dict = {} #the calculated elastic constants go here
    min_eps_dict = {} #here we save the root of the first derivative -> the epsilon at which the energy is minimal 

    for type in deform_types:
        configs = read(f'{directory}{type}.traj@0:', 'r')
        energies = [config.get_potential_energy() for config in configs]
        #[config.get_volume() for config in configs]
        energy_per_volumes = []
        for energy in energies:
            energy_per_volumes.append(energy/volume)

        ax1.plot(epsilons, energy_per_volumes, 'x', label=f'{type}')

        popt, pcov = curve_fit(fit_func, epsilons, energy_per_volumes) 

        ax1.plot(epsilons, [fit_func(x_i, *popt) for x_i in epsilons], '--', label = f'fit {type}')

        derivative_fit_func_0 = lambda eps: derivative_fit_func(eps, *popt[1:]) 

        #find where derivative is 0 -> under what eps energy is minimal
        min_eps = root_scalar(derivative_fit_func_0, x0=0).root
        min_eps_dict[type] = min_eps
        second_derivative_fit_func_0 = lambda eps: second_derivative_fit_func(eps, *popt[2:])
        ax2.plot(epsilons, [derivative_fit_func_0(x_i) for x_i in epsilons], '-')
        ax3.plot(epsilons, [second_derivative_fit_func_0(x_i) for x_i in epsilons], '-')
        elastic_const_dict[type] = second_derivative_fit_func_0(min_eps)
    elastic_const_dict['biax'] =  (elastic_const_dict['biax'] - 2*elastic_const_dict['uniax'])/2
    ax1.legend()
    plt.show()
    fig1.savefig(f'{directory}/strain_epot')
    fig2.savefig(f'{directory}/strain_epot_derivatives')
    #now divide these by lattice const**3
    print(min_eps_dict)
    print(elastic_const_dict)




        

    

