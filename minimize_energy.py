import numpy as np
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.io import Trajectory, read
from ase.eos import EquationOfState

L_sc = 2.56
L_fcc = 3.62
L_bcc = 2.96

sc_cell = SimpleCubic('Cu', latticeconstant=L_sc ,size=(2, 2, 2))
fcc_cell = FaceCenteredCubic('Cu', latticeconstant=L_fcc, size=(2, 2, 2))
bcc_cell = BodyCenteredCubic('Cu', latticeconstant=L_bcc, size=(2,2,2)) 

print(sc_cell.cell[0,0])

traj_sc = Trajectory('sc.traj', 'w')

eps = 0.01
for L in L_sc * np.linspace(1 - eps, 1 + eps, 10):
    sc_cell = SimpleCubic('Cu', latticeconstant=L ,size=(2, 2, 2))
    sc_cell.calc = EMT()
    epot = sc_cell.get_potential_energy()
    print(epot)
    traj_sc.write(sc_cell)

configs = read('sc.traj@0:')
'''
in the next list comp I forgot () after .get_potential_energy()
... no wonder nothing worked
'''
energies = [point.get_potential_energy() for point in configs] 

with open("config_out", 'w') as file:
        for entry in energies:
            file.write(entry)



