
from cc3d import CompuCellSetup




from OP_Cancer_2DSteppables import InitializeAttributesSteppable

CompuCellSetup.register_steppable(steppable=InitializeAttributesSteppable(frequency=1))        




from OP_Cancer_2DSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from OP_Cancer_2DSteppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from OP_Cancer_2DSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from OP_Cancer_2DSteppables import UpdateTumourCellsSteppable

CompuCellSetup.register_steppable(steppable=UpdateTumourCellsSteppable(frequency=1))




from OP_Cancer_2DSteppables import UpdateCAFsSteppable

CompuCellSetup.register_steppable(steppable=UpdateCAFsSteppable(frequency=1))




from OP_Cancer_2DSteppables import UpdateCD8TCellsSteppable

CompuCellSetup.register_steppable(steppable=UpdateCD8TCellsSteppable(frequency=1))




from OP_Cancer_2DSteppables import CD8TCellsMoveSteppable

CompuCellSetup.register_steppable(steppable=CD8TCellsMoveSteppable(frequency=1))




from OP_Cancer_2DSteppables import PlotsSteppable

CompuCellSetup.register_steppable(steppable=PlotsSteppable(frequency=1))



CompuCellSetup.run()
