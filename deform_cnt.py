
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit, root_scalar
from ase.calculators.emt import EMT
from ase.build import nanotube
from ase.io import Trajectory, read, write, gromacs
from ase.optimize import FIRE
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase import units
import os
from pathlib import Path
#from tqdm import tqdm

plot_dir = 'data_CNT_sim/nve_plots/'
directory = 'data_CNT_sim/'
data_directory = 'data_CNT_sim/nve_deform/'
cnt = nanotube(5, 5, symbol='C', length=3, bond=1.42)
cnt.calc = EMT()
#adjust cell and pbc (for visualization)
cnt.cell[:] = cnt.cell[:] + np.array([[10, 0, 0],[0, 10, 0],[0, 0, 0]])
cnt.set_pbc([0,0,1])

start_traj = Trajectory(f'{directory}cnt_initial.traj', 'w')
start_traj.write(cnt)
gromacs.write_gromacs(f'{directory}cnt_initial.gro', cnt) #.gro to lammps.gro so ovito can show it

observing_traj = Trajectory(f'{directory}cnt.traj', 'w')

def deform_in_z(cnt, eps):
    cnt_array = cnt.cell[:]
    cnt_array_new = cnt_array + np.dot(cnt_array_original, (np.array([[0,0,0],[0,0,0],[0,0,eps]])))
    cnt.set_cell(cnt_array_new, scale_atoms=True)

def printenergy(atoms=cnt):
     epot = atoms.get_potential_energy()
     ekin = atoms.get_kinetic_energy()     
     temp = atoms.get_temperature()
     print(f'epot: {epot}, ekin: {ekin}, temperature: {temp}')

def run_simulation(loop_number, step_number, timestep, eps, start_deform=0):
    dyn.attach(printenergy)
    MaxwellBoltzmannDistribution(cnt, temperature_K=0, force_temp=True)
    Stationary(cnt)
    dyn.run(500) #run before simulation (w/o recording data)
    #eqiliberate befor running
    for i in range(loop_number):
    
        loop_time = i*step_number
        for step in range(start_deform): #for if we want to safe data before deforming
                dyn.run(1)
                if (loop_time+step)%1 == 0:
                    data['time'].append((loop_time+step)*timestep)
                    data['temperature'].append(cnt.get_temperature())
                    data['epot'].append(cnt.get_potential_energy())
                    data['ekin'].append(cnt.get_kinetic_energy())
                    data['etot'].append(cnt.get_total_energy())
        deform_in_z(cnt, eps)
        #dyn.attach(printenergy)
        for step in range(start_deform, step_number): 
            #MaxwellBoltzmannDistribution(cnt, temp_K=0)
            #Stationary(cnt)
            dyn.run(1)
            if (loop_time+step)%1 == 0:
                data['time'].append((loop_time+step)*timestep)
                data['temperature'].append(cnt.get_temperature())
                data['epot'].append(cnt.get_potential_energy())
                data['ekin'].append(cnt.get_kinetic_energy())
                data['etot'].append(cnt.get_total_energy())
                if i%1 == 0:
                    gromacs.write_gromacs(f'{data_directory}cnt_{i}.gro', cnt) # dump to lammps.gro so ovito can show it
            
        print(f'Rep {i+1} of {loop_number}')
        if i%10 == 0:
            gromacs.write_gromacs(f'{data_directory}cnt_{i}.gro', cnt) # dump to lammps.gro so ovito can show it
        with open(f'{data_directory}data.json', 'w') as fp:
            json.dump(data, fp)

if __name__ == '__main__':
    for f in Path(f'{data_directory}').glob('*.dump' and '*.gro'): #clean up output path for dump files
        try:
              f.unlink()
        except:
             print('no files in data dir')

    
    strain_rate = 1*10**14 #s^-1 
    #strain_rate = eps/(timestep * step_number*10**-15)
    start_deform = 0
    timestep = 0.1 #in fs
    step_number = 60
    loop_number = 1
    eps = strain_rate * timestep * (step_number-start_deform)*10**-15
    #eps = 0.001
    total_strain = eps*loop_number
    dyn = VelocityVerlet(cnt, timestep=timestep*units.fs)
    data = {'time':[],'temperature':[], 'epot':[],'ekin':[], 'etot': []}
    cnt_array_original = cnt.cell[:]
    #run_simulation(loop_number, step_number, timestep, eps, start_deform) 
    with open(f'{data_directory}data.json', 'r') as fp:
        data = json.load(fp)
    
    fig, ax = plt.subplots()
    ax.set_title('strain rate {a:e} $s^-1$, total strain {b:4f} '.format(a=strain_rate, b=total_strain))
    ax.set_xlabel('time in fs')
    ax.set_ylabel('Epot, Ekin, Etot in eV')
    ax.plot(data['time'], data['temperature'],marker='.', markersize=0.1,  label = 'temperature')
    ax.plot(data['time'], data['ekin'],marker='.', markersize=0.1, label = 'ekin')
    ax.plot(data['time'], data['epot'],marker='.', markersize=0.1, label = 'epot')
    ax.plot(data['time'], data['etot'], marker='.', markersize=0.1, label = 'etot')
    ax.legend()
    plt.show()
    fig.savefig('{plot_dir}long_sr{sr:2.2e}_etot{etot:2.2e}_steps{steps}_loops{loops}.png'.format(plot_dir=plot_dir, sr=strain_rate,
                                                                                    etot=total_strain, steps=step_number,
                                                                                    loops=loop_number))
    print('done')


        





    

