
def configure_simulation():
    
    import csv
    
    # Open file
    
    #parameters_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\parameters.csv" # On DRAC
    parameters_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\parameters.csv" # On local

    # Read parameters

    with open(parameters_file, newline='') as f:
        
        reader = csv.reader(f)
        
        for i, line in enumerate(reader):
            
            if i == 1:
                ifn_diffusion = float(line[1])
            elif i == 2:
                ifn_decay = float(line[1])
            elif i == 3:
                ifn_initial_conc = float(line[1])
            elif i == 4:
                tgf_diffusion = float(line[1])
            elif i == 5:
                tgf_decay = float(line[1])
            elif i == 6:
                tgf_initial_conc = float(line[1])
            elif i == 7:
                collagen_initial_conc = float(line[1])
            
    f.close()
    
    
    from cc3d.core.XMLUtils import ElementCC3D
    
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"5","Version":"4.7.0"})
    
    MetadataElmnt=CompuCell3DElmnt.ElementCC3D("Metadata")
    MetadataElmnt.ElementCC3D("NumberOfProcessors",{},"8")
    MetadataElmnt.ElementCC3D("DebugOutputFrequency",{},"10")
    MetadataElmnt.ElementCC3D("MCSConversionFactor",{"DisplayName":"No conversion","Units":"-"},"1.0")
    MetadataElmnt.ElementCC3D("VoxelConversionFactor",{"DisplayName":"No conversion","Units":"-"},"1.0")
    
    
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    PottsElmnt.ElementCC3D("Dimensions",{"x":"100","y":"100","z":"1"}) # 800 x 1420
    PottsElmnt.ElementCC3D("Steps",{},"500")
    PottsElmnt.ElementCC3D("Temperature",{},"3")
    PottsElmnt.ElementCC3D("NeighborOrder",{},"3")
    
    
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"1","TypeName":"Tumour"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"2","TypeName":"CAF"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"3","TypeName":"myCAF"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"4","TypeName":"CD8T"})
    
    
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Volume"})
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"ExternalPotential"})
    
    
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Tumour"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Medium","Type2":"CAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Medium","Type2":"myCAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Medium","Type2":"CD8T"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Tumour","Type2":"Tumour"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Tumour","Type2":"CAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Tumour","Type2":"myCAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"Tumour","Type2":"CD8T"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"CAF","Type2":"CAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"CAF","Type2":"myCAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"CAF","Type2":"CD8T"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"CD8T","Type2":"CD8T"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"CD8T","Type2":"myCAF"},"10.0")
    PluginElmnt_1.ElementCC3D("Energy",{"Type1":"myCAF","Type2":"myCAF"},"10.0")
    PluginElmnt_1.ElementCC3D("NeighborOrder",{},"1")
    
    
    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"FocalPointPlasticity"})
    
    
    ParametersElmnt=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"Tumour","Type2":"Tumour"})
    
    
    ParametersElmnt.ElementCC3D("Lambda",{},"10")
    ParametersElmnt.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_1=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"Tumour","Type2":"CAF"})
    
    
    ParametersElmnt_1.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_1.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_1.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_1.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_1.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_2=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"Tumour","Type2":"myCAF"})
    ParametersElmnt_2.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_2.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_2.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_2.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_2.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_3=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"Tumour","Type2":"CD8T"})
    ParametersElmnt_3.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_3.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_3.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_3.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_3.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_4=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"CAF","Type2":"CAF"})
    ParametersElmnt_4.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_4.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_4.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_4.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_4.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_5=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"CAF","Type2":"myCAF"})
    ParametersElmnt_5.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_5.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_5.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_5.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_5.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_6=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"CAF","Type2":"CD8T"})
    ParametersElmnt_6.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_6.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_6.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_6.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_6.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_7=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"myCAF","Type2":"myCAF"})
    ParametersElmnt_7.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_7.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_7.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_7.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_7.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_8=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"myCAF","Type2":"CD8T"})
    ParametersElmnt_8.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_8.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_8.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_8.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_8.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    ParametersElmnt_9=PluginElmnt_2.ElementCC3D("Parameters",{"Type1":"CD8T","Type2":"CD8T"})
    ParametersElmnt_9.ElementCC3D("Lambda",{},"10")
    ParametersElmnt_9.ElementCC3D("ActivationEnergy",{},"-50")
    ParametersElmnt_9.ElementCC3D("TargetDistance",{},"7")
    ParametersElmnt_9.ElementCC3D("MaxDistance",{},"20")
    ParametersElmnt_9.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":"1"},"1")
    
    
    PluginElmnt_2.ElementCC3D("NeighborOrder",{},"1")
    
    
    PluginElmnt_3=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"ConnectivityGlobal"})
    PluginElmnt_3.ElementCC3D("Penalty",{"Type":"Tumour"},"1000000")
    PluginElmnt_3.ElementCC3D("Penalty",{"Type":"CAF"},"1000000")
    PluginElmnt_3.ElementCC3D("Penalty",{"Type":"myCAF"},"1000000")
    PluginElmnt_3.ElementCC3D("Penalty",{"Type":"CD8T"},"1000000")
    
    
    #SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"DiffusionSolverFE_OpenCL"})
    SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"DiffusionSolverFE"})
    
    
    DiffusionFieldElmnt=SteppableElmnt.ElementCC3D("DiffusionField",{"Name":"IFN_gamma"})
    DiffusionDataElmnt=DiffusionFieldElmnt.ElementCC3D("DiffusionData")
    DiffusionDataElmnt.ElementCC3D("FieldName",{},"IFN_gamma")
    DiffusionDataElmnt.ElementCC3D("GlobalDiffusionConstant",{},ifn_diffusion)
    DiffusionDataElmnt.ElementCC3D("GlobalDecayConstant",{},ifn_decay)
    DiffusionDataElmnt.ElementCC3D("InitialConcentrationExpression",{},ifn_initial_conc)
    DiffusionFieldElmnt.ElementCC3D("SecretionData")
    
    
    BoundaryConditionsElmnt=DiffusionFieldElmnt.ElementCC3D("BoundaryConditions")
    
    
    PlaneElmnt=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"X"})
    PlaneElmnt.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    PlaneElmnt_1=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"Y"})
    PlaneElmnt_1.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_1.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    
    
    DiffusionFieldElmnt_1=SteppableElmnt.ElementCC3D("DiffusionField",{"Name":"TGF_beta"})
    DiffusionDataElmnt_1=DiffusionFieldElmnt_1.ElementCC3D("DiffusionData")
    DiffusionDataElmnt_1.ElementCC3D("FieldName",{},"TGF_beta")
    DiffusionDataElmnt_1.ElementCC3D("GlobalDiffusionConstant",{},tgf_diffusion)
    DiffusionDataElmnt_1.ElementCC3D("GlobalDecayConstant",{},tgf_decay)
    DiffusionDataElmnt_1.ElementCC3D("InitialConcentrationExpression",{},tgf_initial_conc)
    DiffusionFieldElmnt_1.ElementCC3D("SecretionData")
    
    
    BoundaryConditionsElmnt_1=DiffusionFieldElmnt_1.ElementCC3D("BoundaryConditions")
    
    
    PlaneElmnt_2=BoundaryConditionsElmnt_1.ElementCC3D("Plane",{"Axis":"X"})
    PlaneElmnt_2.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_2.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    PlaneElmnt_3=BoundaryConditionsElmnt_1.ElementCC3D("Plane",{"Axis":"Y"})
    PlaneElmnt_3.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_3.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    
    
    DiffusionFieldElmnt_2=SteppableElmnt.ElementCC3D("DiffusionField",{"Name":"Collagen"})
    
    
    DiffusionDataElmnt_2=DiffusionFieldElmnt_2.ElementCC3D("DiffusionData")
    DiffusionDataElmnt_2.ElementCC3D("FieldName",{},"Collagen")
    DiffusionDataElmnt_2.ElementCC3D("GlobalDiffusionConstant",{},"0")
    DiffusionDataElmnt_2.ElementCC3D("GlobalDecayConstant",{},"0")
    DiffusionDataElmnt_2.ElementCC3D("InitialConcentrationExpression",{},collagen_initial_conc)
    
    
    DiffusionFieldElmnt_2.ElementCC3D("SecretionData")
    
    
    BoundaryConditionsElmnt_2=DiffusionFieldElmnt_2.ElementCC3D("BoundaryConditions")
    
    
    PlaneElmnt_4=BoundaryConditionsElmnt_2.ElementCC3D("Plane",{"Axis":"X"})
    PlaneElmnt_4.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_4.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})
    
    
    PlaneElmnt_5=BoundaryConditionsElmnt_2.ElementCC3D("Plane",{"Axis":"Y"})
    PlaneElmnt_5.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":"0.0"})
    PlaneElmnt_5.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":"0.0"})


    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)    



from cc3d import CompuCellSetup

configure_simulation()            


##############################
## CLASSES FOR INITIALIZING ##
##############################

from OP_Cancer_2DSteppables import HelperFunctionsSteppable
CompuCellSetup.register_steppable(steppable=HelperFunctionsSteppable(frequency=1))          

from OP_Cancer_2DSteppables import InitializeCellPositionSteppable
CompuCellSetup.register_steppable(steppable=InitializeCellPositionSteppable(frequency=1))


#################################
## CLASSES FOR OUTPUTTING DATA ##
#################################

from OP_Cancer_2DSteppables import OutputCSVSteppable
CompuCellSetup.register_steppable(steppable=OutputCSVSteppable(frequency=1))


##################################
## CLASSES FOR BASIC MECHANISMS ##
##################################   


from OP_Cancer_2DSteppables import GrowthSteppable
CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))

from OP_Cancer_2DSteppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))


##########################################
## CLASSES FOR UPDATING CELL PROPERTIES ##
##########################################


from OP_Cancer_2DSteppables import UpdateTumourCellsSteppable
CompuCellSetup.register_steppable(steppable=UpdateTumourCellsSteppable(frequency=1))

from OP_Cancer_2DSteppables import UpdateCAFsSteppable
CompuCellSetup.register_steppable(steppable=UpdateCAFsSteppable(frequency=1))

from OP_Cancer_2DSteppables import UpdateCD8TCellsSteppable
CompuCellSetup.register_steppable(steppable=UpdateCD8TCellsSteppable(frequency=1))


###############################
## CLASSES FOR CELL MOVEMENT ##
###############################


from OP_Cancer_2DSteppables import CD8TCellsMoveSteppable
CompuCellSetup.register_steppable(steppable=CD8TCellsMoveSteppable(frequency=1))

from OP_Cancer_2DSteppables import TumourCellsMoveSteppable
CompuCellSetup.register_steppable(steppable=TumourCellsMoveSteppable(frequency=1))

from OP_Cancer_2DSteppables import CAFsMoveSteppable
CompuCellSetup.register_steppable(steppable=CAFsMoveSteppable(frequency=1))


########################################     
## CLASSES FOR CALIBRATING PARAMETERS ##
########################################

from OP_Cancer_2DSteppables import CellMigrationTrackerSteppable
CompuCellSetup.register_steppable(steppable=CellMigrationTrackerSteppable(frequency=1))

from OP_Cancer_2DSteppables import CD8TKillAttemptsTrackerSteppable
CompuCellSetup.register_steppable(steppable=CD8TKillAttemptsTrackerSteppable(frequency=1))






CompuCellSetup.run()