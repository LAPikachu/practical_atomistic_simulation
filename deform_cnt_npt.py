
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.build import nanotube
from ase.io import Trajectory, read, write, lammpsdata
from ase.optimize import FIRE
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
import os
from pathlib import Path
#from tqdm import tqdm

directory = 'data_CNT_sim/'
data_directory = 'data_CNT_sim/nvt_deform/'
cnt = nanotube(5, 5, symbol='C', length=1)
cnt.calc = EMT()

start_traj = Trajectory(f'{directory}cnt_initial.traj', 'w')
start_traj.write(cnt)
lammpsdata.write_lammps_data(f'{directory}cnt_initial.dump', cnt) # dump to lammps dump so ovito can show it

observing_traj = Trajectory(f'{directory}cnt.traj', 'w')

def deform_in_z(cnt, eps):
    cnt_array = cnt.cell[:]
    cnt_array_new = cnt_array + np.dot(cnt_array_original, (np.array([[0,0,0],[0,0,0],[0,0,eps]])))
    cnt.set_cell(cnt_array_new, scale_atoms=True)

def printenergy(atoms=cnt):
     epot = atoms.get_potential_energy()
     ekin = atoms.get_kinetic_energy()     
     temp = atoms.get_temperature()
     print(f'epot: {epot}, ekin: {ekin}, temperaure: {temp}')

def run_simulation(loop_number, step_number, timestep, eps):
    for i in range(loop_number):
            deform_in_z(cnt, eps)
            #MaxwellBoltzmannDistribution(cnt, temperature_K=300)
            #dyn.attach(printenergy)
            dyn.run(steps=step_number)
            data['time'].append(i*step_number*timestep)
            data['temperatur'].append(cnt.get_temperature())
            data['epot'].append(cnt.get_potential_energy())
            data['ekin'].append(cnt.get_kinetic_energy())
            data['etot'].append(cnt.get_total_energy())
            
            print(f'Rep {i+1} of {loop_number}')
            if i%10 == 0:
                lammpsdata.write_lammps_data(f'{data_directory}cnt_{i}.dump', cnt) # dump to lammps dump so ovito can show it
                with open(f'{data_directory}data.json', 'w') as fp:
                    json.dump(data, fp)

if __name__ == '__main__':
    for f in Path(f'{data_directory}').glob('*.dump'): #clean up output path for dump files
        try:
              f.unlink()
        except:
             print('no files in data dir')

    
    strain_rate = 8 * 10**8
    #strain_rate = eps/(timestep * step_number*10**-12)
    timestep = 0.1 
    step_number = 10
    loop_number = 500
    eps = strain_rate * timestep * step_number*10**-12
    #eps = 0.001
    total_strain = eps*loop_number 
    dyn = Langevin(cnt, timestep=timestep*units.fs, temperature_K= 300, friction=0.01/units.fs)
    data = {'time':[],'temperatur':[], 'epot':[],'ekin':[], 'etot': []}
    cnt_array_original = cnt.cell[:]
    #dyn.attach(printenergy)
    #dyn.run(600)
    #run_simulation(loop_number, step_number, timestep, eps) 
    with open(f'{data_directory}data.json', 'r') as fp:
        data = json.load(fp)
    
    fig, ax = plt.subplots()
    ax.set_title('strain rate {a:e} $s^-1$, total strain {b:4f} '.format(a=strain_rate, b=total_strain))
    ax.set_xlabel('time in fs')
    ax.set_ylabel('Epot, Ekin, Etot in eV, T in K')
    #ax.plot(data['time'], data['temperatur'],'s', label = 'temperatur')
    ax.plot(data['time'], data['ekin'],'o', label = 'ekin')
    ax.plot(data['time'], data['epot'],'x', label = 'epot')
    ax.plot(data['time'], data['etot'], '.', label = 'etot')
    ax.legend()
    plt.show()
    print('done')


        





    

