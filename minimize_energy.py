import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read

#dictionary of initial values for lattice constants (all in Angstrom)
L_init_values = {'sc' : 2.56,
                 'fcc' : 3.62,
                 'bcc' : 2.96}

#traj_sc = Trajectory('sc.traj', 'w')

def create_cubic_lattice(lattice_type, L):
    switch_dict = {'sc' : SimpleCubic('Cu', latticeconstant=L, size=(2,2,2)),
                   'bcc' : BodyCenteredCubic('Cu', latticeconstant=L, size=(2,2,2)),
                   'fcc' : FaceCenteredCubic('Cu', latticeconstant=L, size=(2,2,2))}  
    return switch_dict[lattice_type]

def calculate_latticeconst_epot_traj(lattice_type, L_init):
    traj_file = Trajectory(f'{lattice_type}.traj', 'w')
    lattice_constants_list = []
    for L in np.linspace(L_init - 0.5 , L_init + 1, 100):
        lattice_constants_list.append(L)
        data = create_cubic_lattice(lattice_type, L)
        data.calc = EMT()
        data.get_potential_energy()
        traj_file.write(data)
    return lattice_constants_list

if __name__ =='__main__':
    
    dict = {}
    # loop through crystal structures and lattice constants
    for lattice_type in L_init_values.keys():
        L_init = L_init_values[lattice_type]
        lattice_constants_list = calculate_latticeconst_epot_traj(lattice_type, L_init)
        configs = read(f'{lattice_type}.traj@0:', 'r') #@0: means go though all the saved configurations
        epots_list = [data.get_potential_energy() for data in configs]
        dict[lattice_type] = [lattice_constants_list, epots_list]
    #plot the data
    fig, ax = plt.subplots()
    ax.set_xlabel('Lattice constant in Angstrom')
    ax.set_ylabel('Potential energy in eV')
    for lattice_type in dict.keys():
        x,y = dict[lattice_type]
        ax.plot(x, y, 'o' , label = f'{lattice_type} lattice')
    ax.legend()
    plt.show()
    #fit it to a 2nd order polynomial

    def W(L, p_0, p_1, p_2):
        return p_0 + p_1* L + p_2 * L**2
    
    for lattice_type in L_init_values.keys():
        x,y = dict[lattice_type]
        popt, pcov = curve_fit(W, x, y) 

        ax.plot(x, [W(x_i, *popt) for x_i in x],
                '--', label = f'fit {lattice_type} lattice')


