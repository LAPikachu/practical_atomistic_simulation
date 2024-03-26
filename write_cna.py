from ovito.io import import_file, export_file
from ovito.modifiers import SelectTypeModifier, CommonNeighborAnalysisModifier

pipeline = import_file("data_CNT_sim/nve_deform/cnt_0.gro")

# Insert a CNA modifier to determine the structural type of each atom:
pipeline.modifiers.append(CommonNeighborAnalysisModifier())

# Apply the SelectTypeModifier to select all atoms of FCC and HCP and other type:
pipeline.modifiers.append(SelectTypeModifier(
    operate_on = "particles",
    property = "Structure Type",
    types = { CommonNeighborAnalysisModifier.Type.FCC,
              CommonNeighborAnalysisModifier.Type.HCP,
              CommonNeighborAnalysisModifier.Type.OTHER }
))

# The SelectTypeModifier reports the number of selected elements as an attribute:
data = pipeline.compute()
export_file(data, "output.ss", "imd", columns = 
            ["Particle Identifier", "Particle Type", "Structure type", "Position.X", "Position.Y", "Position.Z"])
print("Number of FCC/HCP atoms: %i" % data.attributes['SelectType.num_selected'])
