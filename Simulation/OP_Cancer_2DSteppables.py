# Imports

from cc3d.core.PySteppables import *

import csv
import math
import numpy as np
import random
from scipy.spatial import KDTree

# SIMULATION VARIABLES

scale = 9 # Micrometers squared per pixel

# GLOBAL VARIABLES

tumour_vol = int(246.3/scale) # Tumour cell area: 246.3 micrometers squared
caf_vol = int(4596/scale) # CAF cell area: 4596.2 micrometers squared
cd8t_vol = int(73.6/scale) # CD8 T cell volume: 73.6 micrometers squared

tumour_lambda_vol = 10
caf_lambda_vol = 10
cd8t_lambda_vol = 50 

tumour_growth = 0.5
caf_growth = 0.05

tumour_apoptosis_probability = 0.00001
caf_apoptosis_probability = 0.00001
cd8t_apoptosis_probability = 0.00001

cd8t_ifn_secretion_rate = 0.1
tumour_tgf_secretion_rate = 0.1
collagen_secretion_rate = 0.1
caf_tgf_secretion_rate = 0.1
caf_ifn_secretion_rate = 0.1

tumour_speed = 1000
caf_speed = 1000
default_cd8t_speed = 1000

caf_activation_threshold_1 = 10000
caf_activation_threshold_2 = 5000
caf_activation_threshold_3 = 1000
caf_activation_prob_1 = 0.01
caf_activation_prob_2 = 0.01
caf_activation_prob_3 = 0.01
caf_activation_prob_4 = 0.01

ifn_pdl1_threshold = 0.05
exhaustion_threshold = 17.5
collagen_threshold = 0.5

# Seed cells based on csv file
cell_position_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\v4_patient28_adata_css.csv"

# Seed cells randomly
total_cell_count = 40


#'''
tumour_proportion = 0.61355
caf_proportion = 0.33880
cd8t_proportion = 0.04765
'''

tumour_proportion = 0.499999
caf_proportion = 0.000001
cd8t_proportion = 0.5
'''

tumour_cd274_proportion = 0.07119
caf_cd274_proportion = 0.11358
cd8t_cd274_proportion = 0.08830

class HelperFunctionsSteppable(SteppableBasePy):
    def in_radius(self, x, y, z, field_type, volume, field):
        dims = field_type.getDim()
        
        radius = math.ceil(math.sqrt(volume / 3.14159))
        lattice_sites = []
        
        # Ensure all cells release the same amount of a given field regardless of size
        
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                if dx**2 + dy**2 <= radius**2:
                    nx, ny, nz = x + dx, y + dy, z
                    if 0 <= nx < dims.x and 0 <= ny < dims.y and 0 <= nz < dims.z:
                        lattice_sites.append((nx, ny, nz))
                        
        for nx, ny, nz in lattice_sites:
            if hasattr(field, 'targetVolume'): # Check if "field" is a cell
                field_type[nx, ny, nz] = field
            else:
                field_type[nx, ny, nz] += (field / len(lattice_sites))


class InitializeCellPositionSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
                
    def start(self):
        
        # Seed cells based on csv file
        '''
        with open(cell_position_file, newline='') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                
                # Set position                
                x = int(float(row["x_aligned"]))
                y = int(float(row["y_aligned"]))
                z = 0
                
                # Set attributes by cell type
                
                cell_type_str = row["leiden_r06"]
                if cell_type_str == "CAF":
                    cell = self.newCell(self.CAF)
                    cell.targetVolume = caf_vol
                    cell.lambdaVolume = caf_lambda_vol
                elif cell_type_str == "Tumour epithelial" or cell_type_str == "Tumour epithelial (proliferative)":
                    cell = self.newCell(self.TUMOUR)
                    cell.targetVolume = tumour_vol
                    cell.lambdaVolume = tumour_lambda_vol
                elif cell_type_str == "CD8 T cell":
                    cell = self.newCell(self.CD8T)
                    cell.targetVolume = cd8t_vol
                    cell.lambdaVolume = cd8t_lambda_vol
                    cell.dict["exhaustion_counter"] = 0
                    cell.dict["speed"] = default_cd8t_speed
                                    
                # Set gene expression
                if float(row["CD274"]) == 0:
                    cell.dict["CD274?"] = False
                else:
                    cell.dict["CD274?"] = True
                
                # Spawn cell
                self.helper_func.in_radius(x, y, z, self.cellField, cell.targetVolume, cell)
        '''
        
        # Seed cells based on random chance
        
        dims = self.cellField.getDim()
               
        # CAF cells
        for i in range(0, int(total_cell_count * caf_proportion)):
            x = random.random() * dims.x
            y = random.random() * dims.y
            z = random.random() * dims.z
            
            cell = self.newCell(self.CAF)
            cell.targetVolume = caf_vol
            cell.lambdaVolume = caf_lambda_vol
            
            if random.random() <= caf_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.in_radius(x, y, z, self.cellField, cell.targetVolume, cell)
        
        # Tumour cells
        for i in range(0, int(total_cell_count * tumour_proportion)):
            x = random.random() * dims.x
            y = random.random() * dims.y
            z = random.random() * dims.z
            
            cell = self.newCell(self.TUMOUR)
            cell.targetVolume = tumour_vol
            cell.lambdaVolume = tumour_lambda_vol
            
            if random.random() <= tumour_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.in_radius(x, y, z, self.cellField, cell.targetVolume, cell)
            
        # CD8T cells
        for i in range(0, int(total_cell_count * cd8t_proportion)):
            x = random.random() * dims.x
            y = random.random() * dims.y
            z = random.random() * dims.z
            
            cell = self.newCell(self.CD8T)
            cell.targetVolume = cd8t_vol
            cell.lambdaVolume = cd8t_lambda_vol
            cell.dict["exhaustion_counter"] = 0
            cell.dict["speed"] = default_cd8t_speed
            
            if random.random() <= cd8t_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.in_radius(x, y, z, self.cellField, cell.targetVolume, cell)
            
        #'''
        
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
    
        for cell in self.cell_list:
            if cell.type == self.TUMOUR: # Tumour cells
                cell.targetVolume += tumour_growth
            elif cell.type == self.CAF or cell.type == self.ACTIVATED_CAF:
                cell.targetVolume += caf_growth    

        
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


    def update_attributes(self):
        
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        
        if self.parent_cell.type == self.TUMOUR or self.parent_cell.type == self.CAF:
            self.child_cell.type = self.parent_cell.type
        elif self.parent_cell.type == self.ACTIVATED_CAF:
            self.child_cell.type = self.CAF
    
        # reinitialise dict attributes based on cell type
        if self.child_cell.type == self.CD8T:  # CD8T
            self.child_cell.dict["exhaustion_counter"] = 0
            self.child_cell.dict["cd8t-speed"] = default_cd8t_speed
            self.child_cell.dict["CD274?"] = False
        else:  # Tumour, CAF, or activated CAF
            self.child_cell.dict["CD274?"] = False

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        

class UpdateTumourCellsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
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
            cd8t = None
            
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(tumour):
                if neighbor:
                    if neighbor.type == self.CD8T: #
                        cd8t = neighbor

                        self.helper_func.in_radius(cd8t.xCOM, cd8t.yCOM, cd8t.zCOM, self.field.IFN_gamma,
                            tumour.targetVolume, cd8t_ifn_secretion_rate)
                        self.shared_steppable_vars["total_ifn_gamma"] += cd8t_ifn_secretion_rate
                        
                        # If both are not expressing CD274 and CD8T cell is not at exhaustion threshold, CD8 T cell successfully kills tumour cell:
                        if (cd8t.dict["CD274?"] == False or tumour.dict["CD274?"] == False) and cd8t.dict["exhaustion_counter"] < exhaustion_threshold:
                            cd8t.dict["exhaustion_counter"] += 1
                            cells_to_delete.append(tumour)
                            continue
                        # Else, immune escape occurs
                        else:
                            cd8t.dict["exhaustion_counter"] = exhaustion_threshold
                    
            # Check 2    
            if self.field.IFN_gamma[int(tumour.xCOM), int(tumour.yCOM), int(tumour.zCOM)] > ifn_pdl1_threshold:
                tumour.dict["CD274?"] = True
            
            self.helper_func.in_radius(tumour.xCOM, tumour.yCOM, tumour.zCOM, self.field.TGF_beta,
                tumour.targetVolume, tumour_tgf_secretion_rate)
            self.shared_steppable_vars["total_tgf_beta"] += tumour_tgf_secretion_rate
    
        for tumour in cells_to_delete:
            self.delete_cell(tumour)
            self.shared_steppable_vars["dead_tumour_count"] += 1


class UpdateCAFsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
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
                        cd8t.dict["exhaustion_threshold"] = exhaustion_threshold
        
        # Check 3
        if caf.dict["CD274?"] == False and self.field.IFN_gamma[caf.xCOM, caf.yCOM, caf.zCOM] > 0:
            caf.dict["CD274?"] = True
        
    def step(self, mcs):
        
        cells_to_delete = []
        
        for caf in self.cell_list_by_type(self.CAF): # Inactivated CAF
            
            # Check 2 & 3
            self.all_cafs_checks(caf, cells_to_delete)
                
            # Check 4: Calculate probability of CAF activation
            if self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM] >= caf_activation_threshold_1:
                if random.random() <= caf_activation_prob_1:
                    caf.type = self.ACTIVATED_CAF
            elif self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM] >= caf_activation_threshold_2:
                if random.random() <= caf_activation_prob_2:
                    caf.type = self.ACTIVATED_CAF
            elif self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM] >= caf_activation_threshold_3:
                if random.random() <= caf_activation_prob_3:
                    caf.type = self.ACTIVATED_CAF
            else:
                if random.random() <= caf_activation_prob_4:
                    caf.type = self.ACTIVATED_CAF
                    caf.type = self.ACTIVATED_CAF
                
            # Action
            self.helper_func.in_radius(caf.xCOM, caf.yCOM, caf.zCOM, self.field.TGF_beta,
                caf.targetVolume, caf_tgf_secretion_rate)
            self.helper_func.in_radius(caf.xCOM, caf.yCOM, caf.zCOM, self.field.IFN_gamma,
                caf.targetVolume, caf_ifn_secretion_rate)
        
        for a_caf in self.cell_list_by_type(self.ACTIVATED_CAF):
            
            # "Check" 1
            self.helper_func.in_radius(a_caf.xCOM, a_caf.yCOM, a_caf.zCOM, self.field.Collagen,
                a_caf.targetVolume, collagen_secretion_rate)
            
            # Check 2 & 3
            self.all_cafs_checks(a_caf, cells_to_delete)
            
            # Action
            self.helper_func.in_radius(a_caf.xCOM, a_caf.yCOM, a_caf.zCOM, self.field.TGF_beta,
                a_caf.targetVolume, caf_tgf_secretion_rate)
            self.shared_steppable_vars["total_tgf_beta"] += caf_tgf_secretion_rate
            
        for caf in cells_to_delete:
            self.delete_cell(caf)
            self.shared_steppable_vars["dead_caf_count"] += 1
    

class UpdateCD8TCellsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.helper_func = HelperFunctionsSteppable()
        
    def start(self):
        self.shared_steppable_vars["dead_cd8t_count"] = 0
        self.shared_steppable_vars["exhausted_cd8t"] = 0
               
    def step(self, mcs):
        
        cells_to_delete = []
        
        # Apoptosis rate
        for cd8t in self.cell_list_by_type(self.CD8T):
            if random.random() <= cd8t_apoptosis_probability:
                cells_to_delete.append(cd8t)
                continue
        
        for cd8t in self.cell_list_by_type(self.CD8T):
        
            # Check 1 done in UpdateTumourCellsSteppable and UpdateCAFsSteppable
            # Check 2 done in UpdateTumourCellsSteppable
                
            # Check 3 & 4
                        
            if self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM] > collagen_threshold:
                cd8t.dict["speed"] = 0
            else:
                cd8t.dict["speed"] = default_cd8t_speed
                
        for cd8t in cells_to_delete:
            self.delete_cell(cd8t)
            self.shared_steppable_vars["dead_cd8t_count"] +=1
            

class TumourCellsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
    
    def start(self):
        
        for tumour in self.cell_list_by_type(self.TUMOUR):
            tumour.lambdaVecX = tumour_speed * uniform(-0.5,0.5)
            tumour.lambdaVecY = tumour_speed * uniform(-0.5,0.5)
    
    def step(self, mcs):
        
        for tumour in self.cell_list_by_type(self.TUMOUR):
            tumour.lambdaVecX = tumour_speed * uniform(-0.5,0.5)
            tumour.lambdaVecY = tumour_speed * uniform(-0.5,0.5)

            
class CAFsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
    
    def start(self):
        
        for caf in (self.cell_list_by_type(self.CAF) + self.cell_list_by_type(self.Activated_CAF)):
            caf.lambdaVecX = caf_speed * uniform(-0.5,0.5)
            caf.lambdaVecY = caf_speed * uniform(-0.5,0.5)
    
    def step(self, mcs):
        
        for caf in (self.cell_list_by_type(self.CAF) + self.cell_list_by_type(self.Activated_CAF)):
            caf.lambdaVecX = caf_speed * uniform(-0.5,0.5)
            caf.lambdaVecY = caf_speed * uniform(-0.5,0.5)
        

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
            
            norm = (dx**2 + dy**2)**0.5
            
            if norm > 0:
                cd8t.lambdaVecX = dx/norm * cd8t.dict["speed"]
                cd8t.lambdaVecY = dy/norm * cd8t.dict["speed"]
             
            
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
        tumour_expressing_cd274 = sum(1 for tumour in self.cell_list_by_type(self.TUMOUR)
            if tumour.dict["CD274?"] == True)
        all_caf_expressing_cd274 = sum(1 for caf in self.cell_list_by_type(self.CAF) if caf.dict["CD274?"] == True) + \
            sum(1 for a_caf in self.cell_list_by_type(self.ACTIVATED_CAF) if a_caf.dict["CD274?"] == True)
        cd8t_expressing_cd274 = sum(1 for cd8t in self.cell_list_by_type(self.CD8T) 
            if cd8t.dict["CD274?"] == True)
        
        self.plot_gene_expression.add_data_point("Tumour With CD274",
            mcs, tumour_expressing_cd274)
        self.plot_gene_expression.add_data_point("CAF With CD274",
            mcs, all_caf_expressing_cd274)
        self.plot_gene_expression.add_data_point("CD8T With CD274",
            mcs, cd8t_expressing_cd274)
            