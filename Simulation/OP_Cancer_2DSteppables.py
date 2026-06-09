# Imports

from cc3d.core.PySteppables import *

import csv
import math
import numpy as np
import random
from scipy.spatial import KDTree

# Global variables

# Model parameters
tumour_vol = 9
cd8t_vol = 9
caf_vol = 9

tumour_growth = 0.1 # Get from Sameena
caf_growth = 0.05

tumour_apoptosis_probability = 0.001
caf_apoptosis_probability = 0.001

kill_ifn_secretion_rate = 0.1
tumour_tgf_secretion_rate = 0.1
collagen_secretion_rate = 0.1
caf_tgf_secretion_rate = 0.1
caf_ifn_secretion_rate = 0.1

default_cd8t_speed = 100

caf_activation_threshold = 0.5
exhaustion_threshold = 35
caf_radius = 10
cd8t_slowdown_rate = 3
collagen_threshold = 0.5

tumour_pdl1_expression = 0.8
caf_pdl1_expression = 0.8
cd8t_pdcd1_expression = 0.8

cell_position_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\v4_patient28_adata_css.csv"

class InitializeCellPositionSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def seed_cell_at(self, cell, x, y, z):
        
        radius = int(math.sqrt(cell.targetVolume / 3.14))
        dims = self.cellField.getDim()
        
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                nx, ny, nz = x + dx, y + dy, z
                if 0 <= nx < dims.x and 0 <= ny < dims.y and 0 <= nz < dims.z:
                    self.cellField[nx, ny, nz] = cell
            
        
    def start(self):
        
        with open(cell_position_file, newline='') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                
                # Set position                
                x = int(float(row["x_aligned"]))
                y = int(float(row["y_aligned"]))
                z = 0
                
                # Set cell type
                cell_type_str = row["leiden_r06"]
                if cell_type_str == "Tumour epithelial" or cell_type_str == "Tumour epithelial (proliferative)":
                    cell_type = getattr(self, "TUMOUR")
                elif cell_type_str == "CD8 T cell":
                    cell_type = getattr(self, "CD8T")
                elif cell_type_str == "CAF":
                    cell_type = getattr(self, "CAF")
                
                cell = self.newCell(cell_type)
                
                if cell.type == self.TUMOUR:
                    cell.targetVolume = tumour_vol
                    cell.lambdaVolume = 10
                elif cell.type == self.CAF or cell.type == self.ACTIVATED_CAF:
                    cell.targetVolume = caf_vol
                    cell.lambdaVolume = 10
                elif cell.type == self.CD8T:
                    cell.targetVolume = cd8t_vol
                    cell.lambdaVolume = 15
                    
                # Set gene expression
                if float(row["CD274"]) == 0:
                    cell.dict["CD274?"] = False
                else:
                    cell.dict["CD274?"] = True
                    
                # Set other attributes
                if cell.type == self.CD8T:
                    cell.dict["exhausted?"] = False
                    cell.dict["exhaustion_time"] = 0
                    cell.dict["speed"] = default_cd8t_speed
                
                # Spawn cell
                self.seed_cell_at(cell, x, y, z)


class InitializeAttributesSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        '''
        for cell in self.cell_list:
            if cell.type == self.TUMOUR:
                
                if random.random() <= tumour_pdl1_expression:
                    cell.dict["CD274?"] = True
                else:
                    cell.dict["CD274?"] = False
                
            elif cell.type == self.CAF:
                
                if random.random() <= caf_pdl1_expression:
                    cell.dict["CD274?"] = True
                else:
                    cell.dict["CD274?"] = False
                    
            elif cell.type == self.ACTIVATED_CAF:
                
                if random.random() <= caf_pdl1_expression:
                    cell.dict["CD274?"] = True
                else:
                    cell.dict["CD274?"] = False

            elif cell.type == self.CD8T:
                cell.dict["exhausted?"] = False
                cell.dict["exhaustion_time"] = 0
                cell.dict["speed"] = default_cd8t_speed
                
                if random.random() <= cd8t_pdcd1_expression:
                    cell.dict["CD274?"] = True
                else:
                    cell.dict["CD274?"] = False
        '''        
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
    
        for cell in self.cell_list:
            if cell.type == self.TUMOUR: # Tumour cells
                cell.targetVolume += tumour_growth
            elif cell.type == self.CAF or cell.type == self.ACTIVATED_CAF:
                cell.targetVolume += caf_growth

        # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
            # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # you can use here any fcn of concentrationAtCOM
            # cell.targetVolume += 0.01 * concentrationAtCOM       

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list:
            if cell.type == self.TUMOUR:
                if cell.volume > 2*tumour_vol:
                    cells_to_divide.append(cell)
            elif cell.type == self.CAF or cell.type == self.ACTIVATED_CAF:
                if cell.volume > 2*caf_vol:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        
        if self.parent_cell.type == self.TUMOUR or self.parent_cell.type == self.CAF:
            self.child_cell.type = self.parent_cell.type
        elif self.parent_cell.type == self.ACTIVATED_CAF:
            self.child_cell.type = self.CAF
    
        # reinitialise dict attributes based on cell type
        if self.child_cell.type == self.CD8T:  # CD8T
            self.child_cell.dict["exhausted"] = False
            self.child_cell.dict["exhaustion_time"] = 0
            self.child_cell.dict["cd8t-speed"] = 1
            self.child_cell.dict["CD274?"] = False
        else:  # Tumour, CAF, or activated CAF
            self.child_cell.dict["CD274?"] = False

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        

class UpdateTumourCellsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        self.shared_steppable_vars["dead_tumour_count"] = 0
        self.shared_steppable_vars["total_ifn_gamma"] = 0
        self.shared_steppable_vars["total_tgf_beta"] = 0
        
    def step(self, mcs):
        
        cells_to_delete = []
        
        for i, tumour in enumerate(self.cell_list_by_type(self.TUMOUR)):
            
            # Check for apoptosis
            if random.random() <= tumour_apoptosis_probability:
                cells_to_delete.append(tumour)
                continue
            
            # Check 1
            immune_escaped = False
            cd8t = None
            
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(tumour):
                if neighbor:
                    if neighbor.type == self.CD8T: #
                        cd8t = neighbor
                        if cd8t.dict["CD274?"] == True and tumour.dict["CD274?"] == True:
                            cd8t.dict["exhausted?"] = True
                            immune_escaped = True      
            
            if cd8t != None:
                self.field.IFN_gamma[tumour.xCOM, tumour.yCOM, tumour.zCOM] += kill_ifn_secretion_rate
                self.shared_steppable_vars["total_ifn_gamma"] += caf_ifn_secretion_rate
                if not immune_escaped:
                    cells_to_delete.append(tumour)
                else:
                    
                    # Check 2
                    
                    if self.field.IFN_gamma[int(tumour.xCOM), int(tumour.yCOM), int(tumour.zCOM)] > 0:
                        tumour.dict["CD274?"] = True
                        
                    self.field.TGF_beta[tumour.xCOM, tumour.yCOM, tumour.zCOM] += tumour_tgf_secretion_rate
                    self.shared_steppable_vars["total_tgf_beta"] += tumour_tgf_secretion_rate
    
        for tumour in cells_to_delete:
            self.delete_cell(tumour)
            self.shared_steppable_vars["dead_tumour_count"] += 1


class UpdateCAFsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        self.shared_steppable_vars["dead_caf_count"] = 0
        
    def all_cafs_checks(self, caf, cells_to_delete):
        
        # Check for death
        if random.random() <= caf_apoptosis_probability:
            cells_to_delete.append(caf)
            return
        
        # Check 2
        cd8t = None
            
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(caf):
            if neighbor:
                if neighbor.type == self.CD8T:
                    cd8t = neighbor
                    if cd8t.dict["CD274?"] == True and caf.dict["CD274?"] == True:
                        cd8t.dict["exhausted"] = True
        
        # Check 3
        if caf.dict["CD274?"] == False and self.field.IFN_gamma[caf.xCOM, caf.yCOM, caf.zCOM] > 0:
            caf.dict["CD274?"] = True
        
    def step(self, mcs):
        
        cells_to_delete = []
        
        for caf in self.cell_list_by_type(self.CAF):
            
            # Check 2 & 3
            self.all_cafs_checks(caf, cells_to_delete)
                
            # Check 4
            if self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM] > caf_activation_threshold:
                caf.type = self.ACTIVATED_CAF
                
            # Action
            self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM] += caf_tgf_secretion_rate
            self.field.IFN_gamma[caf.xCOM, caf.yCOM, caf.zCOM] += caf_ifn_secretion_rate
        
        for a_caf in self.cell_list_by_type(self.ACTIVATED_CAF):
            
            # "Check" 1
            self.field.Collagen[a_caf.xCOM, a_caf.yCOM, a_caf.zCOM] += collagen_secretion_rate
            
            # Check 2 & 3
            self.all_cafs_checks(a_caf, cells_to_delete)
            
            # Action
            self.field.TGF_beta[a_caf.xCOM, a_caf.yCOM, a_caf.zCOM] += caf_tgf_secretion_rate
            self.shared_steppable_vars["total_tgf_beta"] += caf_tgf_secretion_rate
            
        for caf in cells_to_delete:
            self.delete_cell(caf)
            self.shared_steppable_vars["dead_caf_count"] += 1
    

class UpdateCD8TCellsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
    def start(self):
        self.shared_steppable_vars["dead_cd8t_count"] = 0
        self.shared_steppable_vars["exhausted_cd8t"] = 0
               
    def step(self, mcs):
        
        cells_to_delete = []
        
        all_CAFs = list(self.cell_list_by_type(self.CAF)) + list(self.cell_list_by_type(self.ACTIVATED_CAF))
        CAF_positions = [(caf.xCOM, caf.yCOM, caf.zCOM) for caf in all_CAFs ]
        
        if len(CAF_positions) == 0:
            return
        
        CAF_tree = KDTree(CAF_positions)
        
        for cd8t in self.cell_list_by_type(self.CD8T):
        
            # Check 1 done in UpdateTumourCellsSteppable and UpdateCAFsSteppable
            
            # Check 2
            if cd8t.dict["exhausted?"] == True:
                cd8t.dict["exhaustion_time"] += 1            
            if cd8t.dict["exhaustion_time"] > exhaustion_threshold:
                cells_to_delete.append(cd8t)
                continue
                
            # Check 3 & 4
            
            neighbours = CAF_tree.query_ball_point((cd8t.xCOM, cd8t.yCOM, cd8t.zCOM), caf_radius)
            
            if len(neighbours) > 0:
                cd8t.dict["speed"] = max([0, cd8t.dict["speed"] - (len(neighbours) * cd8t_slowdown_rate)])
            elif self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM] > collagen_threshold:
                cd8t.dict["speed"] = 0
            else:
                cd8t.dict["speed"] = default_cd8t_speed
                
        for cd8t in cells_to_delete:
            self.delete_cell(cd8t)
            self.shared_steppable_vars["dead_cd8t_count"] +=1
            


class CD8TCellsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
    
    def step(self, mcs):
        
        tumour_cells = list(self.cell_list_by_type(self.TUMOUR))
        
        if len(tumour_cells) == 0:
            return
        
        tumour_positions = [(tumour.xCOM, tumour.yCOM, tumour.zCOM) for tumour in tumour_cells ]
        tumour_tree = KDTree(tumour_positions)
         
        for cd8t in self.cell_list_by_type(self.CD8T):
            
            distance, index = tumour_tree.query((cd8t.xCOM, cd8t.yCOM, cd8t.zCOM))
            nearest_tumour = tumour_cells[index]
            
            dx = cd8t.xCOM - nearest_tumour.xCOM
            dy = cd8t.yCOM - nearest_tumour.yCOM
            
            cd8t.lambdaVecX = dx * cd8t.dict["speed"]
            cd8t.lambdaVecY = dy * cd8t.dict["speed"]
             
            
class PlotsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):        
        
        # Plot live cell count
        self.plot_live_count = self.add_new_plot_window(title='Live Cell Count over Time',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count',
                                                 x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False,
                                                 config_options={'legend':True})
        
        self.plot_live_count.add_plot("Live Tumour Count", style='Lines', color='coral', size=5)       
        self.plot_live_count.add_plot("Live CAF Count", style='Lines', color='yellowgreen', size=5)
        self.plot_live_count.add_plot("Live Activated CAF Count", style='Lines', color='gold', size=5)
        self.plot_live_count.add_plot("Live CD8T Count", style='Lines', color='cornflowerblue', size=5)
       
        # Plot dead cell count
        self.plot_dead_count = self.add_new_plot_window(title='Dead Cell Count over Time',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count',
                                                 x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False,
                                                 config_options={'legend':True})
        
        self.plot_dead_count.add_plot("Dead Tumour Count", style='Lines', color='coral', size=5)       
        self.plot_dead_count.add_plot("Dead CD8T Count", style='Lines', color='cornflowerblue', size=5)
        self.plot_dead_count.add_plot("Dead CAF Count", style='Lines', color='yellowgreen', size=5)
        
        # Plot total diffusing agent
        self.plot_diffusing_agent = self.add_new_plot_window(title='Total Diffusing Agent over Time',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Total Diffusing Agent',
                                                 x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False,
                                                 config_options={'legend':True})
        
        self.plot_diffusing_agent.add_plot("Total IFN-Gamma", style='Lines', color='lightsteelblue', size=5)       
        self.plot_diffusing_agent.add_plot("Total TGF-Beta", style='Lines', color='rosybrown', size=5)
        
        # Plot gene expression
        self.plot_gene_expression = self.add_new_plot_window(title='Cell Count With Gene Expression',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Cell Count',
                                                 x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False,
                                                 config_options={'legend':True})
        
        self.plot_gene_expression.add_plot("Tumour With CD274",
            style='Lines', color='lightsteelblue', size=5)       
        self.plot_gene_expression.add_plot("CAF With CD274",
            style='Lines', color='rosybrown', size=5)
        self.plot_gene_expression.add_plot("CD8T With CD274",
            style='Lines', color='burlywood', size=5)
        
    def step(self,mcs):
  
        tumour_count = 0; caf_count = 0; activated_caf_count = 0; cd8t_count = 0;
        
        # Plot live cell count
        for cell in self.cell_list:
                if cell.type == self. TUMOUR:
                    tumour_count += 1
                elif cell.type == self.CAF:
                    caf_count += 1
                elif cell.type == self.ACTIVATED_CAF:
                    activated_caf_count += 1
                elif cell.type == self.CD8T:
                    cd8t_count += 1

        self.plot_live_count.add_data_point("Live Tumour Count", mcs, tumour_count)
        self.plot_live_count.add_data_point("Live CAF Count", mcs, caf_count)
        self.plot_live_count.add_data_point("Live Activated CAF Count", mcs, activated_caf_count)
        self.plot_live_count.add_data_point("Live CD8T Count", mcs, cd8t_count)
        
        # Plot dead cell count
        self.plot_dead_count.add_data_point("Dead Tumour Count",
            mcs, self.shared_steppable_vars["dead_tumour_count"])
        self.plot_dead_count.add_data_point("Dead CD8T Count",
            mcs, self.shared_steppable_vars["dead_cd8t_count"])
        self.plot_dead_count.add_data_point("Dead CAF Count",
            mcs, self.shared_steppable_vars["dead_caf_count"])
          
        # Plot total diffusing agent
        self.plot_diffusing_agent.add_data_point("Total IFN-Gamma", mcs, self.shared_steppable_vars["total_ifn_gamma"])
        self.plot_diffusing_agent.add_data_point("Total TGF-Beta", mcs, self.shared_steppable_vars["total_tgf_beta"])
        
        # Plot gene expression
        tumour_expressing_pdl1 = sum(1 for tumour in self.cell_list_by_type(self.TUMOUR)
            if tumour.dict["CD274?"] == True)
        all_caf_expressing_pdl1 = sum(1 for caf in self.cell_list_by_type(self.CAF) if caf.dict["CD274?"] == True) + \
            sum(1 for a_caf in self.cell_list_by_type(self.ACTIVATED_CAF) if a_caf.dict["CD274?"] == True)
        cd8t_expressing_pdcd1 = sum(1 for cd8t in self.cell_list_by_type(self.CD8T) 
            if cd8t.dict["CD274?"] == True)
        
        self.plot_gene_expression.add_data_point("Tumour With CD274",
            mcs, tumour_expressing_pdl1)
        self.plot_gene_expression.add_data_point("CAF With CD274",
            mcs, all_caf_expressing_pdl1)
        self.plot_gene_expression.add_data_point("CD8T With CD274",
            mcs, cd8t_expressing_pdcd1)
            