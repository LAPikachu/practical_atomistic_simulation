import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.visualize import view
from matplotlib.widgets import Slider

radius = 0.3
fig, axs = plt.subplots(1, 3, figsize=(15, 1))
fig.subplots_adjust(left=0.25, bottom=0.25)
sc_cell = SimpleCubic('Cu', latticeconstant=1 ,size=(2, 2, 2))
fcc_cell = FaceCenteredCubic('Cu', latticeconstant=1 ,size=(2, 2, 2))
bcc_cell = BodyCenteredCubic('Cu', latticeconstant=1, size=(2,2,2)) 

def make_plot(radius):
    for i in range(len(axs)): #clear every plot seprately
        axs[i].cla()
    plot_atoms(sc_cell, ax=axs[0], radii=radius)
    plot_atoms(fcc_cell, ax=axs[1], radii=radius)
    plot_atoms(bcc_cell, ax=axs[2], radii=radius)
    axs[0].set_title('Simple Centered Cubic Lattice')
    axs[1].set_title('Face-Centered Cubic Lattice')
    axs[2].set_title('Body-Centered Cubic Lattice')

make_plot(radius)
axslider = fig.add_axes([0.1, 0.25, 0.025, 0.63])

radius_slider = Slider(
    ax=axslider,
    label="radius",
    valmin=0.0,
    valmax=1.0,
    valinit=radius,
    orientation="vertical"
)

def update(val):
    make_plot(radius_slider.val)
    fig.canvas.draw_idle()

radius_slider.on_changed(update)
plt.show()
fig.savefig("crystal_lattices.png")


#view(fcc_cell)
