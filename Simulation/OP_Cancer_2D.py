
from cc3d import CompuCellSetup


###################################################
## CLASSES FOR SIMULATING BIOPHYSICAL MECHANISMS ##
###################################################

from OP_Cancer_2DSteppables import HelperFunctionsSteppable
CompuCellSetup.register_steppable(steppable=HelperFunctionsSteppable(frequency=1))          

from OP_Cancer_2DSteppables import InitializeCellPositionSteppable
CompuCellSetup.register_steppable(steppable=InitializeCellPositionSteppable(frequency=1))          

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

from OP_Cancer_2DSteppables import TumourCellsMoveSteppable
CompuCellSetup.register_steppable(steppable=TumourCellsMoveSteppable(frequency=1))


########################################     
## CLASSES FOR CALIBRATING PARAMETERS ##
########################################

from OP_Cancer_2DSteppables import CellSpeedTrackerSteppable
CompuCellSetup.register_steppable(steppable=CellSpeedTrackerSteppable(frequency=1))

from OP_Cancer_2DSteppables import CD8TKillAttemptsTrackerSteppable
CompuCellSetup.register_steppable(steppable=CD8TKillAttemptsTrackerSteppable(frequency=1))


#################################
## CLASSES FOR OUTPUTTING DATA ##
#################################

from OP_Cancer_2DSteppables import OutputCSVSteppable
CompuCellSetup.register_steppable(steppable=OutputCSVSteppable(frequency=1))

from OP_Cancer_2DSteppables import PlotsSteppable
CompuCellSetup.register_steppable(steppable=PlotsSteppable(frequency=1))



CompuCellSetup.run()
