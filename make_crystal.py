import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
slab = FaceCenteredCubic('Cu', size=(2, 2, 2))
fig, ax = plt.subplots()
plot_atoms(slab, ax, radii=0.3)
fig.savefig("ase_slab.png")
view(slab)