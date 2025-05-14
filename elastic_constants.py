import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
#from ase.calculators.emt import EMT
from asap3 import EMT, EMTRasmussenParameters, EMT2013
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read
from ase.optimize import FIRE

#print(cell.cell[:]) # cell.cell only gives a 1-dim array type output :/ ?

strain_high = 0.1

strain_low = -0.1

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
    atoms_positions_new = np.dot(atoms_positions,(np.identity(3)+strain_tensor))
    cell_array_new = np.dot(cell_array,(np.identity(3)+strain_tensor))
    cell.cell[:] = cell_array_new
    cell.positions = atoms_positions_new
    return cell

def fit_func(eps, p_0, p_1, p_2, p_3):
        return p_0 + p_1*eps + p_2*eps**2 + p_3*eps**3
    
def derivative_fit_func(eps, p_1, p_2, p_3):
    return p_1 + 2*p_2*eps + 3*p_3*eps**2

def second_derivative_fit_func(eps,p_2,p_3):
    return 2*p_2 + 6*p_3*eps

if __name__ == '__main__':
    #relax initial structur
    with open("data_min_energy/lattice_consts", "r") as file: #load lattice constant
        lattice_const_dict = json.load(file) 
    lattice_const = lattice_const_dict['W_2nd_ord']['fcc']
   
    cell = FaceCenteredCubic('Cu', latticeconstant=lattice_const, size=(3,3,3))
    cell.calc = EMT()
    # cell.calc = EMT(EMTRasmussenParameters()) in case other parameters need to be used
    #relax cell  
    traj_relaxed_cell = Trajectory(f'{directory}relaxed_cell.log', 'w')
    relax_cell =  FIRE(cell, trajectory=f'{directory}relaxed_cell.traj')
    relax_cell.run(fmax=10**-8)
    traj_relaxed_cell.write(cell)
    traj_relax_read = Trajectory(f'{directory}relaxed_cell.traj', 'r')
    cell_relaxed = traj_relax_read[-1]

    
    #get initial state of relaxed cell to work with
    cell_array = cell_relaxed.cell[:]
    atoms_positions = cell_relaxed.positions
    volume = cell_relaxed.get_volume()

    epsilons = np.linspace(strain_low, strain_high, 60)
    deform_types = ['uniax', 'biax', 'shear']
    #deform the cell 
    for d_type in deform_types:
        cell = read(f'{directory}relaxed_cell.traj@0', 'r')
        #traj = Trajectory(f'{directory}{d_type}.traj', 'r')
        #cell = traj[0]
        tensor = StrainTensor()
        traj = Trajectory(f'{directory}{d_type}.traj', 'w')

        for eps in epsilons:
            tensor.make_tensor(d_type, eps)
            cell = deform_cell(cell_array, tensor.tensor)
            relax_cell.run(fmax=10**-8)
            cell.calc = EMT()
            cell.get_potential_energy()
            traj.write(cell)
    
    fig1, ax1 = plt.subplots()
    fig2, (ax2, ax3) = plt.subplots(1,2)
    ax1.set_xlabel('Strain $\epsilon$')
    ax1.set_ylabel('Strain Energy density eV/$Angstrom^{3}$')
    ax2.set_ylabel('$dW/d\epsilon$ eV/$Angstrom^{3}$')
    ax3.set_ylabel('$d^{2}W/d\epsilon^{2}$ eV')
    ax2.set_xlabel('Strain $\epsilon$')
    ax3.set_xlabel('Strain $\epsilon$')
    elastic_const_dict = {} #the calculated elastic constants go here
    min_eps_dict = {} #here we save the root of the first derivative -> the epsilon at which the energy is minimal 
    
    #fit and find min
    for d_type,mark in zip(deform_types,['x','s','o']):
        configs = read(f'{directory}{d_type}.traj@0:', 'r')
        energies = [config.get_potential_energy() for config in configs]
        #[config.get_volume() for config in configs]
        energy_per_volumes = []
        for energy in energies:
            energy_per_volumes.append(energy/volume)

        ax1.plot(epsilons, energy_per_volumes, mark , label=f'{d_type}')

        popt, pcov = curve_fit(fit_func, epsilons, energy_per_volumes) 

        ax1.plot(epsilons, [fit_func(x_i, *popt) for x_i in epsilons], '--', label = f'fit {d_type}')

        derivative_fit_func_0 = lambda eps: derivative_fit_func(eps, *popt[1:]) 

        #find where derivative is 0 -> under what eps energy is minimal
        min_eps = root_scalar(derivative_fit_func_0, x0=0).root
        min_eps_dict[d_type] = min_eps
        second_derivative_fit_func_0 = lambda eps: second_derivative_fit_func(eps, *popt[2:])
        ax2.plot(epsilons, [derivative_fit_func_0(x_i) for x_i in epsilons], '-')
        ax3.plot(epsilons, [second_derivative_fit_func_0(x_i) for x_i in epsilons], '-')
        elastic_const_dict[d_type] = second_derivative_fit_func_0(min_eps)*1.6021 #1eV/Angstrom = 1.6021 N/m**2
    elastic_const_dict['biax'] =  (elastic_const_dict['biax'] - 2*elastic_const_dict['uniax'])/2 #correction for C_12
    elastic_const_dict['shear'] = elastic_const_dict['shear']/4
    ax1.legend()
    plt.show()
    fig1.savefig(f'{directory}/strain_epot_elow{strain_low}_ehigh{strain_high}.png')
    fig2.savefig(f'{directory}/strain_epot_derivatives_elow{strain_low}_ehigh{strain_high}.png')
    #now divide these by lattice const**3
    with open(f'{directory}elastic_constants.json', 'w') as fp:
        json.dump(elastic_const_dict, fp)
    print(min_eps_dict)
    print(elastic_const_dict)

    print('done')




        

    

