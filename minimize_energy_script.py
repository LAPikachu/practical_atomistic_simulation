# the following lines include all the imports we need

import numpy as np
import matplotlib.pyplot as plt 
import json
import csv
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
#from asap3 import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read
import os
from tqdm import tqdm, trange

'''
List of crystal lattice types
sc: simple cubic
fcc: face centered cubic
bcc: body centered cubic
'''
lattice_types = ['sc', 'fcc', 'bcc']

#dictionary of initial values for lattice constants (all in Angstrom)
#in task 3.1 you estimated values for the lattice constants, insert in the dictionary the below
L_init_values = {'sc' : 1, 
                 'fcc' : 1,
                 'bcc': 1}

def create_cubic_lattice(lattice_type, L):
    switch_dict = {'sc' : SimpleCubic('Cu', latticeconstant=L, size=(1,1,1)),
                   'bcc' : BodyCenteredCubic('Cu', latticeconstant=L, size=(1,1,1)),
                   'fcc' : FaceCenteredCubic('Cu', latticeconstant=L, size=(1,1,1))}  
    return switch_dict[lattice_type]

def calculate_latticeconst_epot_traj(lattice_type, L_init):
    os.makedirs(os.path.dirname('/data_min_energy'), exist_ok=True)
    traj_file = Trajectory(f'data_min_energy/{lattice_type}.traj', 'w')
    lattice_constants_list = []
    L_high = L_init + 0.3 
    L_low = L_init - 0.3 
    N = 300
    iterator_values = np.linspace(L_low, L_high, N)
    print(f'\n\n\nLattice type : {lattice_type}')
    print(f"Calculating {N} data point for lattice constants {L_low:.2f} to {L_high:.2f}")
    for L in tqdm(iterator_values):
        lattice_constants_list.append(L) # append current lattice constant
        data = create_cubic_lattice(lattice_type, L)
        data.calc = EMT()
        data.get_potential_energy()
        traj_file.write(data)
    return lattice_constants_list

def write_data_dict(lattice_type ,lattice_constants_list):
    data_dict = {}
    L_init = L_init_values[lattice_type]
    os.makedirs(os.path.dirname('/data_min_energy'), exist_ok=True)
    configs = read(f'data_min_energy/{lattice_type}.traj@0:', 'r') #@0:  go though all the saved configurations
    epots_list = [data.get_potential_energy() for data in configs] #list potential energies
    data_dict[lattice_type] = [lattice_constants_list, epots_list] #array pairing epot w/ corresponding latticeconst
    return data_dict[lattice_type]

#save data from dicts to permanent file
def save_data_to_json(data_dict, filename):

    with open(filename, 'w') as data:
        json.dump(data_dict, data)

def write_to_csv(lattice_type, strains, energy_densities):
    os.makedirs(os.path.dirname('/data_min_energy'), exist_ok=True)
    with open(f'data_min_energy/{lattice_type}_strain_energydensity.csv', 'w', newline='') as new_csv:
        newarray = csv.writer(new_csv, delimiter=",")
        newarray.writerow(["strain in angstrom", "energy density in eV/Angstrom^3"])
        data =  []
        for strain, energy_density in zip(strains, energy_densities):
            data.append([strain, energy_density])
        newarray.writerows(data)

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
    print("Give initial values for lattice constants in Angstom (e.g. 1.9) \n")
    L_init_values["sc"] = float(input("Simple Cubic: "))
    L_init_values["fcc"] = float(input("Face Centered Cubic: "))
    L_init_values["bcc"] = float(input("Body Centered Cubic: "))
    data_dict = {}
    for lattice_type in lattice_types:
        L_init = L_init_values[lattice_type]
        lattice_constants_list = calculate_latticeconst_epot_traj(lattice_type, L_init)
        data_dict[lattice_type] = write_data_dict(lattice_type, lattice_constants_list)
    os.makedirs(os.path.dirname('/data_min_energy'), exist_ok=True)
    save_data_to_json(data_dict, 'data_min_energy/epot_latticeconst')

    with open('data_min_energy/epot_latticeconst', 'r') as data:
       data_dict = json.load(data)
    print("\n\n\n")
    for lattice_type in lattice_types:
        strains, energy_densities = data_dict[lattice_type]
        write_to_csv(lattice_type, strains, energy_densities)
    #plot the data
        print(f"Strain vs energy density written to csv-files to file 'data_min_energy/{lattice_type}_strain_energydensity.csv'")