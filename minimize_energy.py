import numpy as np
import matplotlib.pyplot as plt 
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read

#Task 3

#List of crystal lattice types
lattice_types = ['sc', 'fcc', 'bcc']

#dictionary of initial values for lattice constants (all in Angstrom)
L_init_values = {'sc' : 2.56,
                 'fcc' : 3.62,
                 'bcc' : 2.96}

def create_cubic_lattice(lattice_type, L):
    switch_dict = {'sc' : SimpleCubic('Cu', latticeconstant=L, size=(1,1,1)),
                   'bcc' : BodyCenteredCubic('Cu', latticeconstant=L, size=(1,1,1)),
                   'fcc' : FaceCenteredCubic('Cu', latticeconstant=L, size=(1,1,1))}  
    return switch_dict[lattice_type]

def calculate_latticeconst_epot_traj(lattice_type, L_init):
    traj_file = Trajectory(f'data_min_energy/{lattice_type}.traj', 'w')
    lattice_constants_list = []
    for L in np.linspace(L_init - 0.2 , L_init + 0.2, 100):
        lattice_constants_list.append(L) # append current lattice constant
        data = create_cubic_lattice(lattice_type, L)
        data.set_calculator(EMT())
        data.get_potential_energy()
        traj_file.write(data)
    return lattice_constants_list

def write_data_dict(lattice_type ,lattice_constants_list):
    data_dict = {}
    L_init = L_init_values[lattice_type]
    calculate_latticeconst_epot_traj(lattice_type, L_init)
    configs = read(f'data_min_energy/{lattice_type}.traj@0:', 'r') #@0:  go though all the saved configurations
    epots_list = [data.get_potential_energy() for data in configs] #list potential energies
    data_dict[lattice_type] = [lattice_constants_list, epots_list] #array pairing epot w/ corresponding latticeconst
    return data_dict[lattice_type]

#save data from dicts to permanent file
def save_data_to_json(data_dict, filename):
    with open(filename, 'w') as data:
        json.dump(data_dict, data)

def fit_L_min_epot(func, dfunc, data_dict):
    fitting_params = {}
    func_name = func.__name__
    for lattice_type in lattice_types:
        x,y = data_dict[lattice_type]#[:60][:60] could adjust over which values to fit
        popt, pcov = curve_fit(func, x, y)
        ax1.plot(x, [func(x_i, *popt) for x_i in x],
                '--', label = f'{func_name} {lattice_type} lattice')
        fitting_params[lattice_type] = popt
    ax1.legend()
    fig1.savefig(f'data_min_energy/{func_name}_epot_lattice_plot')
    roots_dict = {}
    for lattice_type in  lattice_types:
        p_array = fitting_params[lattice_type][1:] # not from 0 since p_0 sits at 0
        L_init = L_init_values[lattice_type]
        dfunc_fit = lambda x: dfunc(x, p_array)
        roots_dict[lattice_type] = root_scalar(dfunc_fit, x0=L_init).root #finds min
    return roots_dict, func_name 

if __name__ =='__main__':

    data_dict = {}
    for lattice_type in lattice_types:
        L_init = L_init_values[lattice_type]
        lattice_constants_list = calculate_latticeconst_epot_traj(lattice_type, L_init)
        data_dict[lattice_type] = write_data_dict(lattice_type, lattice_constants_list)
    save_data_to_json(data_dict, 'data_min_energy/epot_latticeconst')


    with open('data_min_energy/epot_latticeconst', 'r') as data:
       data_dict = json.load(data)
    #plot the data
    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Lattice constant in Angstrom')
    ax1.set_ylabel('Potential energy in eV')
    # for plotting
    for lattice_type, marker_type in zip(lattice_types, ['1','<','x']):
        x,y = data_dict[lattice_type]
        ax1.plot(x, y, marker_type , label = f'{lattice_type} lattice')
    ax1.legend()
    #fit it to a 2nd order polynomial
    #for making the method dynamic we initialize p as a list
    def W_2nd_ord(L, p_0, p_1, p_2):
        return p_0 + p_1* L + p_2 * L**2
    #derivative of function
    #p is not a variable since scip.optimiz.root_scalar should not fit over it
    def dW_2nd(L, p_array):
        return p_array[0] + 2 * p_array[1]* L  
    
    roots_dict, func_name = fit_L_min_epot(W_2nd_ord, dW_2nd, data_dict)
    roots_fit_dict = {}
    roots_fit_dict[func_name] = roots_dict
    '''
    fit is not really optimal
    2nd ord poly not ideal for theses kinds of curvature?
    '''
    
    '''
    def W_3rd_ord(L, p_0, p_1, p_2, p_3):
        return p_0 + p_1 * L + p_2 * L ** 2 + p_3 * L ** 3
    def dW_3rd_ord(L, p_array):
        return p_array[0] + 2 * p_array[1] * L + 3 * p_array[2] * L**2
    
    roots_dict, func_name = fit_L_min_epot(W_3rd_ord, dW_3rd_ord, data_dict)
    '''
    roots_fit_dict[func_name] = roots_dict

    save_data_to_json(roots_fit_dict, 'data_min_energy/lattice_consts')
    print("Lattice constants: \n"
        "Sc {L_sc} \nbcc {L_bcc}"
        "\nfcc {L_fcc} ".format(L_sc=roots_dict['sc'],
                                L_bcc=roots_dict['bcc'], L_fcc=roots_dict['fcc']))