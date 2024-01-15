import numpy as np
import matplotlib.pyplot as plt 
import json
from scipy.optimize import curve_fit, root 
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read

#List of crystal lattice types
lattice_types = ['sc', 'fcc', 'bcc']

#dictionary of initial values for lattice constants (all in Angstrom)
L_init_values = {'sc' : 2.56,
                 'fcc' : 3.62,
                 'bcc' : 2.96}

def create_cubic_lattice(lattice_type, L):
    switch_dict = {'sc' : SimpleCubic('Cu', latticeconstant=L, size=(2,2,2)),
                   'bcc' : BodyCenteredCubic('Cu', latticeconstant=L, size=(2,2,2)),
                   'fcc' : FaceCenteredCubic('Cu', latticeconstant=L, size=(2,2,2))}  
    return switch_dict[lattice_type]

def calculate_latticeconst_epot_traj(lattice_type, L_init):
    traj_file = Trajectory(f'{lattice_type}.traj', 'w')
    lattice_constants_list = []
    for L in np.linspace(L_init - 0.5 , L_init + 1, 100):
        lattice_constants_list.append(L) # append current lattice constant
        data = create_cubic_lattice(lattice_type, L)
        data.calc = EMT()
        data.get_potential_energy()
        traj_file.write(data)
    return lattice_constants_list

def write_data_dict(lattice_type ,lattice_constants_list):
    data_dict = {}
    L_init = L_init_values[lattice_type]
    calculate_latticeconst_epot_traj(lattice_type, L_init)
    configs = read(f'{lattice_type}.traj@0:', 'r') #@0:  go though all the saved configurations
    epots_list = [data.get_potential_energy() for data in configs] #list potential energies
    data_dict[lattice_type] = [lattice_constants_list, epots_list] #array pairing epot w/ corresponding latticeconst
    return data_dict[lattice_type]

#save data from dicts to permanent file
def save_data_to_json(data_dict, filename):
    with open(filename, 'w') as data:
        json.dump(data_dict, data)

if __name__ =='__main__':
    '''
    data_dict = {}
    for lattice_type in lattice_types:
        L_init = L_init_values[lattice_type]
        lattice_constants_list = calculate_latticeconst_epot_traj(lattice_type, L_init)
        data_dict[lattice_type] = write_data_dict(lattice_type, lattice_constants_list)
    save_data_to_json(data_dict, 'epot_latticeconst')
    '''

    with open('epot_latticeconst', 'r') as data:
       data_dict = json.load(data)
    #plot the data
    fig, ax = plt.subplots()
    ax.set_xlabel('Lattice constant in Angstrom')
    ax.set_ylabel('Potential energy in eV')
    # for plotting
    for lattice_type, marker_type in zip(lattice_types, ['1','<','x']):
        x,y = data_dict[lattice_type]
        ax.plot(x, y, marker_type , label = f'{lattice_type} lattice')
    ax.legend()
    #fit it to a 2nd order polynomial

    def W_2nd(L, p_0, p_1, p_2):
        return p_0 + p_1* L + p_2 * L**2

    fitting_params = {}
    for lattice_type in lattice_types:
        x,y = data_dict[lattice_type]#[:60][:60] could adjust over which values to fit
        popt, pcov = curve_fit(W_2nd, x, y) 
        ax.plot(x, [W_2nd(x_i, *popt) for x_i in x],
                '--', label = f'2nd order fit {lattice_type} lattice')
        fitting_params[lattice_type] = popt
    ax.legend()
    fig.savefig('epot_lattice_plot')
    '''
    fit is not really optimal
    2nd ord poly not ideal for theses kinds of curvature?
    '''
    roots_dict = {}
    for lattice_type in  lattice_types:
        p_1, p_2 = fitting_params[lattice_type][1:] # not from 0 since p_0 sits at 0
        def dW_2nd_zero(L):
            return p_1 + 2 * p_2 * L
        L_init = L_init_values[lattice_type]
        roots_dict[lattice_type] = root(dW_2nd_zero, L_init).x[0]
    save_data_to_json(roots_dict, 'lattice_consts')
    print("Lattice constants: \n"
          "Sc {L_sc} \nbcc {L_bcc}"
          "\nfcc {L_fcc} ".format(L_sc=roots_dict['sc'],
                                  L_bcc=roots_dict['bcc'], L_fcc=roots_dict['fcc']))