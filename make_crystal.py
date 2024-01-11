import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.visualize import view

sc_cell = SimpleCubic('Cu', latticeconstant=1 ,size=(2, 2, 2))
fcc_cell = FaceCenteredCubic('Cu', latticeconstant=1 ,size=(2, 2, 2))
bcc_cell = BodyCenteredCubic('Cu', latticeconstant=1, size=(2,2,2)) 

fig, axs = plt.subplots(1, 3, figsize=(10, 5))

plot_atoms(sc_cell, ax=axs[0], radii=radius)
plot_atoms(fcc_cell, ax=axs[1], radii=radius)
plot_atoms(bcc_cell, ax=axs[2], radii=radius)



axs[0].set_title('Simple Cubic Lattice')
axs[1].set_title('Face-Centered Cubic Lattice')
axs[2].set_title('Body-Centered Cubic Lattice')
fig.savefig("crystal_lattices.png")

#view(fcc_cell)
