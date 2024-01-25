
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.build import nanotube
from ase.io import Trajectory, read, write, lammpsdata
from ase.optimize import FIRE
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
#from tqdm import tqdm

directory = 'data_CNT_sim/'
cnt = nanotube(5, 5, symbol='C', length=1)
cnt.calc = EMT()

start_traj = Trajectory(f'{directory}cnt_initial.traj', 'w')
start_traj.write(cnt)
lammpsdata.write_lammps_data(f'{directory}cnt_initial.dump', cnt) # dump to lammps dump so ovito can show it

observing_traj = Trajectory(f'{directory}cnt.traj', 'w')

def deform_in_z(cnt, eps):
    cnt_array = cnt.cell[:]
    cntcnt_array_new = np.dot(cnt_array, (np.eye(3) + np.array([[0,0,0],[0,0,0],[0,0,eps]])))
    cnt.cell = cnt_array_new




if __name__ == '__main__':

    MaxwellBoltzmannDistribution(cnt, temperature_K=0)
    dyn = VelocityVerlet(cnt, timestep=0.1*units.fs)
    eps = 0.10
    data = {'time':[],'temperatur':[], 'epot':[],'ekin':[]}
    
    dyn.run(100)

    for i in range(10):
        data['time'].append(i)
        deform_in_z(cnt, eps)
        print(cnt.cell[:])
        dyn.run(10)
        data['temperatur'].append(cnt.get_temperature())
        data['epot'].append(cnt.get_potential_energy())
        data['ekin'].append(cnt.get_kinetic_energy())
        print(f'Rep {i} of 100')
        lammpsdata.write_lammps_data(f'{directory}cnt.dump', cnt) # dump to lammps dump so ovito can show it

    with open(f'{directory}data.json', 'w') as fp:
        json.dump(data, fp)
    
    with open(f'{directory}data.json', 'r') as fp:
        data = json.load(fp)
    
    fig, ax = plt.subplots()
    ax.plot(data['time'], data['temperatur'], label = 'temperatur')
    ax.plot(data['time'], data['ekin'], label = 'ekin')
    ax.plot(data['time'], data['epot'], label = 'epot')
    ax.legend()
    plt.show()
    


        





    

