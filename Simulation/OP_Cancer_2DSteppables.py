# Imports

from cc3d.core.PySteppables import *

import csv
import math
import numpy as np
import os
import random
from random import uniform
from scipy.spatial import KDTree

######################
## GLOBAL VARIABLES ##
######################

# Open files

#parameters_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\parameters.csv" 
#cell_position_file = r"C:\CompuCell3D\Projects\OP_Cancer_2D\patient28_truncated_normalized_filtered.csv" # On local PC

# On DRAC
parameters_file = r"/home/annied/OP_Cancer_2D/parameters.csv"
cell_position_file = r"/home/annied/OP_Cancer_2D/patient28_edited.csv" 

# Read parameters

with open(parameters_file, newline='') as f:
    
    reader = csv.reader(f)
    
    for i, line in enumerate(reader):
        if i == 8:
            tumour_vol = round(float(line[1]))
        elif i == 9:
            caf_vol = round(float(line[1]))
        elif i == 10:
            cd8t_vol = round(float(line[1]))
        elif i == 11:
            tumour_prolif = float(line[1])
        elif i == 12:
            caf_prolif = float(line[1])
        elif i == 13:
            tumour_apoptosis_prob = float(line[1])
        elif i == 14:
            caf_apoptosis_prob = float(line[1])
        elif i == 15:
            cd8t_apoptosis_prob = float(line[1])
        elif i == 16:
            tumour_migration = float(line[1])
        elif i == 17:
            caf_migration = float(line[1])
        elif i == 18:
            default_cd8t_migration = float(line[1])
        elif i == 19:
            cd8t_ifn_secretion = float(line[1])
        elif i == 20:
            tumour_tgf_secretion = float(line[1])
        elif i == 21:
            collagen_secretion = float(line[1])
        elif i == 22:
            caf_tgf_secretion = float(line[1])
        elif i == 23:
            exhaustion_threshold = float(line[1])
        elif i == 24:
            tumour_ifn_pdl1_threshold = float(line[1])
        elif i == 25:
            caf_ifn_pdl1_threshold = float(line[1])
            
f.close()

# Set established parameters (will not change

tumour_lambda_vol = 50 
caf_lambda_vol = 10
cd8t_lambda_vol = 50 

# Seed cells randomly
total_cell_count = 40

tumour_proportion = 0.92307
caf_proportion = 0
cd8t_proportion = 0.0769

tumour_cd274_proportion = 0.07119
caf_cd274_proportion = 0.11358
cd8t_cd274_proportion = 0.08830


##############################
## CLASSES FOR INITIALIZING ##
##############################
  

class HelperFunctionsSteppable(SteppableBasePy):
    def update_lattice_sites(self, x, y, z, field_type, volume, value):
        '''
        If [field_type] is the cell field and [value] is a cell:
            Change [volume] lattice sites of [field_type] to the cell type of [value]
        
        Otherwise:
            Change [volume] lattice sites of [field_type] by a total of [value] around [x,y,z]
        
        Args:
            x, y, z: the central coordinates around which the lattice sites will be modified
            field_type: the field (cell field, TGF-beta, IFN-gamma, or collagen) whose lattice sites will be modified
            volume: the number of lattice sites that will be modified 
            value: the modification to the lattice site
        
        '''
        
        # Find the coordinates of the [volume] lattice sites in a circle around x,y,z  
        dims = field_type.getDim()
        radius = math.ceil(math.sqrt(volume / 3.14159))
        lattice_sites = []
               
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                if dx**2 + dy**2 <= radius**2:
                    nx, ny, nz = x + dx, y + dy, z
                    if 0 <= nx < dims.x and 0 <= ny < dims.y and 0 <= nz < dims.z:
                        lattice_sites.append((nx, ny, nz))
        
        # Modify the values of the lattice sites
        for nx, ny, nz in lattice_sites:
            # Check if "value" is a cell
            # If so, set the cell field to that cell type
            if hasattr(value, 'targetVolume'):
                field_type[nx, ny, nz] = value
            # Else, increment the value of the field
            else:
                # Ensure all cells release the same amount of a given field regardless of size
                field_type[nx, ny, nz] += (value / len(lattice_sites))


class InitializeCellPositionSteppable(SteppableBasePy):
    
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
        # Create helper object
        self.helper_func = HelperFunctionsSteppable()
                
    def start(self):
        
        '''
        Initializes cell type, position, and CD274 expression based on spatial transcriptomics data.
        Initializes cell size based on parameters from literature.
        '''
        
        dims = self.cellField.getDim()
        
        # Seed cells based on csv file
        #'''
        with open(cell_position_file, newline='') as f:
            reader = csv.DictReader(f)
            
            
            for row in reader:
                
                # Set position                
                x = int(float(row["x_aligned"]))
                y = int(float(row["y_aligned"]))
                z = 0
                
                if x >= dims.x or y >= dims.y or z >= dims.z:
                    continue
                
                # Set size by cell type
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
                    cell.dict["migration"] = default_cd8t_migration
                                    
                # Set gene expression
                if float(row["CD274"]) == 0:
                    cell.dict["CD274?"] = False
                else:
                    cell.dict["CD274?"] = True
                    
                # Track movement   
                cell.dict["position_history"] = [x, y, z]
                
                # Spawn cell
                self.helper_func.update_lattice_sites(x, y, z, self.cellField, cell.targetVolume, cell)
                
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
            
            cell.dict["position_history"] = [x, y, z]
            
            if random.random() <= caf_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.update_lattice_sites(x, y, z, self.cellField, cell.targetVolume, cell)
        
        # Tumour cells
        for i in range(0, int(total_cell_count * tumour_proportion)):
            x = random.random() * dims.x
            y = random.random() * dims.y
            z = random.random() * dims.z
            
            cell = self.newCell(self.TUMOUR)
            cell.targetVolume = tumour_vol
            cell.lambdaVolume = tumour_lambda_vol
            
            cell.dict["position_history"] = [x, y, z]
            
            if random.random() <= tumour_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.update_lattice_sites(x, y, z, self.cellField, cell.targetVolume, cell)
            
        # CD8T cells
        for i in range(0, int(total_cell_count * cd8t_proportion)):
            x = random.random() * dims.x
            y = random.random() * dims.y
            z = random.random() * dims.z
            
            cell = self.newCell(self.CD8T)
            cell.targetVolume = cd8t_vol
            cell.lambdaVolume = cd8t_lambda_vol
            cell.dict["exhaustion_counter"] = 0
            cell.dict["migration"] = default_cd8t_migration
            
            cell.dict["position_history"] = [x, y, z]
            
            if random.random() <= cd8t_cd274_proportion:
                cell.dict["CD274?"] = True
            else:
                cell.dict["CD274?"] = False
            
            self.helper_func.update_lattice_sites(x, y, z, self.cellField, cell.targetVolume, cell)
            
        #'''


##################################
## CLASSES FOR BASIC MECHANISMS ##
##################################  

        
class GrowthSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        '''
        Increment cell size based on growth rate.
        '''
        
        for cell in self.cell_list:
            if cell.type == self.TUMOUR: # Tumour cells
                cell.targetVolume *= tumour_prolif
            elif cell.type == self.CAF or cell.type == self.MYCAF:
                cell.targetVolume *= caf_prolif    

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

    def step(self, mcs):
        '''
        Perform mitosis when cell volume is equal to or greater than 2 times its default size.
        '''

        cells_to_divide=[]
        for cell in self.cell_list:
            if cell.type == self.TUMOUR:
                if cell.volume >= 2*tumour_vol:
                    cells_to_divide.append(cell)
            elif cell.type == self.CAF or cell.type == self.MYCAF:
                if cell.volume >= 2*caf_vol:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)


    def update_attributes(self):
        '''
        Initialize attributes of child cells.
        '''
        
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        
        # Initialize cell type
        if self.parent_cell.type == self.TUMOUR or self.parent_cell.type == self.CAF:
            self.child_cell.type = self.parent_cell.type
        elif self.parent_cell.type == self.MYCAF:
            self.child_cell.type = self.CAF
    
        # Initialize attributes based on cell type
        if self.child_cell.type == self.CD8T:  # CD8T
            self.child_cell.dict["exhaustion_counter"] = 0
            self.child_cell.dict["cd8t-migration"] = default_cd8t_migration
            self.child_cell.dict["CD274?"] = False
        else:  # Tumour, CAF, or myCAF
            self.child_cell.dict["CD274?"] = False


##########################################
## CLASSES FOR UPDATING CELL PROPERTIES ##
##########################################  


class UpdateTumourCellsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
    def start(self):
              
        # Track model outputs
        '''
        Move the tumour cell in a random direction at the tumour migration rate.
        '''
        self.shared_steppable_vars["dead_tumour_count"] = 0
        self.shared_steppable_vars["total_ifn_gamma"] = 0
        self.shared_steppable_vars["total_tgf_beta"] = 0
        
        # Track CD8 T kills        
        self.shared_steppable_vars["cd8t_kill_attempts"] = {}
        total_mcs = self.simulator.getNumSteps()
        for cd8t in self.cell_list_by_type(self.CD8T):
            self.shared_steppable_vars["cd8t_kill_attempts"][cd8t.id] = [0, total_mcs]
            
        
    def step(self, mcs):
        '''
        Update attributes of tumour cells and their effects on other cells.
        '''
        
        cells_to_delete = []
        
        for i, tumour in enumerate(self.cell_list_by_type(self.TUMOUR)):
            
            # Baseline apoptosis
            if random.random() <= tumour_apoptosis_prob:
                cells_to_delete.append(tumour)
                continue
            
            # CHECK 1: if neighbouring CD8 T cell
            cd8t = None
            
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(tumour):
                if neighbor:
                    if neighbor.type == self.CD8T:
                        cd8t = neighbor
                        
                        # CD8 T cell secretes IFN gamma
                        self.helper_func.update_lattice_sites(cd8t.xCOM, cd8t.yCOM, cd8t.zCOM, self.field.IFN_gamma,
                            tumour.targetVolume, cd8t_ifn_secretion)
                        self.shared_steppable_vars["total_ifn_gamma"] += cd8t_ifn_secretion
                        
                        # Track CD8 T kills
                        self.shared_steppable_vars["cd8t_kill_attempts"][cd8t.id][0] += 1
                        
                        # Check if CD8 T cell successfully kills tumour cell:
                        if (cd8t.dict["CD274?"] == False or tumour.dict["CD274?"] == False) and cd8t.dict["exhaustion_counter"] < exhaustion_threshold:
                            cd8t.dict["exhaustion_counter"] += 1
                            cells_to_delete.append(tumour)
                            continue
                        # Else, immune escape and CD8 T exhaustion occurs
                        else:
                            cd8t.dict["exhaustion_counter"] = exhaustion_threshold
                    
            # CHECK 2: if CD274 expression is induced
            if self.field.IFN_gamma[int(tumour.xCOM), int(tumour.yCOM), int(tumour.zCOM)] > tumour_ifn_pdl1_threshold:
                tumour.dict["CD274?"] = True
            
            self.helper_func.update_lattice_sites(tumour.xCOM, tumour.yCOM, tumour.zCOM, self.field.TGF_beta,
                tumour.targetVolume, tumour_tgf_secretion)
            self.shared_steppable_vars["total_tgf_beta"] += tumour_tgf_secretion
        
        # Delete tumour cells marked for apoptosis die
        for tumour in cells_to_delete:
            self.delete_cell(tumour)
            self.shared_steppable_vars["dead_tumour_count"] += 1


class UpdateCAFsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
    def start(self):
        
        # Track model outputs
        self.shared_steppable_vars["dead_caf_count"] = 0
        
    def all_cafs_checks(self, caf, cells_to_delete):
        '''
        Update CAFs and myCAFs.
        '''
        
        # Baseline apoptosis
        if random.random() <= caf_apoptosis_prob:
            cells_to_delete.append(caf)
            return
        
        # CHECK 2: if neighbouring CD8 T cells
        cd8t = None
        
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(caf):
            if neighbor:
                if neighbor.type == self.CD8T:
                    cd8t = neighbor
                    # Check if CD8 T celle exhaustion occurs
                    if cd8t.dict["CD274?"] == True and caf.dict["CD274?"] == True:
                        cd8t.dict["exhaustion_threshold"] = exhaustion_threshold
        
        # CHECK 3: if CD274 expression is induced
        if caf.dict["CD274?"] == False and self.field.IFN_gamma[caf.xCOM, caf.yCOM, caf.zCOM] > caf_ifn_pdl1_threshold:
            caf.dict["CD274?"] = True
        
    def step(self, mcs):
        '''
        Update attributes of CAF cells and their effect on other cells.
        '''
        
        cells_to_delete = []
        
        for caf in self.cell_list_by_type(self.CAF):
            
            # CHECK 2 & 3
            self.all_cafs_checks(caf, cells_to_delete)
            
            
            # CHECK 4: if myCAF phenotype induced
            tgf = self.field.TGF_beta[caf.xCOM, caf.yCOM, caf.zCOM]
            prob = 1 - (math.exp((math.log(1-(-0.0029*(tgf**2) + 0.0562*tgf + 0.1084)))/86400))
            
            if random.random() <= prob:                    
                caf.type = self.MYCAF
                
            # Secrete TGF-beta
            self.helper_func.update_lattice_sites(caf.xCOM, caf.yCOM, caf.zCOM, self.field.TGF_beta,
                caf.targetVolume, caf_tgf_secretion)
        
        for mycaf in self.cell_list_by_type(self.MYCAF):
            
            # myCAF secrete collagen
            self.helper_func.update_lattice_sites(mycaf.xCOM, mycaf.yCOM, mycaf.zCOM, self.field.Collagen,
                mycaf.targetVolume, collagen_secretion)
            
            # CHECK 2 & 3
            self.all_cafs_checks(mycaf, cells_to_delete)
            
            # Secrete TGF-beta
            self.helper_func.update_lattice_sites(mycaf.xCOM, mycaf.yCOM, mycaf.zCOM, self.field.TGF_beta,
                mycaf.targetVolume, caf_tgf_secretion)
            self.shared_steppable_vars["total_tgf_beta"] += caf_tgf_secretion
        
        # Delete CAFs marked for apoptosis die
        for caf in cells_to_delete:
            self.delete_cell(caf)
            self.shared_steppable_vars["dead_caf_count"] += 1
        

class UpdateCD8TCellsSteppable(SteppableBasePy):
    
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.helper_func = HelperFunctionsSteppable()
        
    def start(self):
        
        # Track model outputs
        self.shared_steppable_vars["dead_cd8t_count"] = 0
        self.shared_steppable_vars["exhausted_cd8t"] = 0
               
    def step(self, mcs):
        '''
        Update attributes of CD8 T cells and their effects on other cells.
        '''
        
        cells_to_delete = []
        
        # Apoptosis rate
        for cd8t in self.cell_list_by_type(self.CD8T):
            if random.random() <= cd8t_apoptosis_prob:
                cells_to_delete.append(cd8t)
                continue
        
        for cd8t in self.cell_list_by_type(self.CD8T):
        
            # Check 1 done in UpdateTumourCellsSteppable and UpdateCAFsSteppable
            # Check 2 done in UpdateTumourCellsSteppable
                
            # CD8 T cell migration is affected by collagen density           
            collagen = int(self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM])
            
            cd8t.dict["migration"] = default_cd8t_migration
            
                        
            collagen = self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM]
            cd8t.dict["migration"] = -2585 * collagen + 780
            
        # Delete CD8 T cells marked for apoptosis 
        for cd8t in cells_to_delete:
            self.shared_steppable_vars["dead_cd8t_count"] +=1
            
            # Track CD8 T kills
            self.shared_steppable_vars["cd8t_kill_attempts"][cd8t.id][1] = mcs
            
            self.delete_cell(cd8t)
            
            
###############################
## CLASSES FOR CELL MOVEMENT ##
###############################


class CD8TCellsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def compute_shift(self, cell, distance, end):
        '''
        Compute how far CD8 T cells can move based on possible obstructions in their path.
        
        Arguments:
            cell: CD8 T cell that will move.
            distance: the farthest distance the CD8 T cell can move.
            end: the farthest coordinates the CD8 T cell can move to.
            
        Returns:
            The direction vector of the CD8 T cell, indicating how many lattice sites in the direction of each axis (x,y,z) the cell will move.
        '''
                
        dims = self.cellField.getDim()
        
        x0 = cell.xCOM; y0 = cell.yCOM; z0 = cell.zCOM
        x1, y1, z1 = end
        
        dest_x = x0; dest_y = y0; dest_z = z0
        
        for i in range (1, int(round(distance)) + 1):
            t = i / distance
            x = int(round(x0 + (x1 - x0) * t))
            y = int(round(y0 + (y1 - y0) * t))
            z = int(round(z0 + (z1 - z0) * t))
        
            if not (0 <= x < dims.x and 0 <= y < dims.y and 0 <= z < dims.z):
                break
                
            occupant = self.cell_field[x, y, z]
            
            if occupant is None or occupant.id == cell.id:
                dest_x = x; dest_y = y; dest_z = z
            else:
                break
                
        return(int(round(dest_x - x0)), int(round(dest_y - y0)), int(round(dest_z - z0)))
        
    
    def step(self, mcs):
        '''
        Find the tumour cell nearest to each CD8 T cell and move the CD8 T cell in that direction at the appropriate migration rate.
        '''
                       
        tumour_cells = list(self.cell_list_by_type(self.TUMOUR))
        
        if len(tumour_cells) == 0:
            return
        
        tumour_positions = [(tumour.xCOM, tumour.yCOM, tumour.zCOM) for tumour in tumour_cells ]
        tumour_tree = KDTree(tumour_positions)
         
        for cd8t in self.cell_list_by_type(self.CD8T):
            
            distance, index = tumour_tree.query((cd8t.xCOM, cd8t.yCOM, cd8t.zCOM))
            nearest_tumour = tumour_cells[index]
                                
            if distance > 0:     
                                              
                # Set new pixel
                if distance <= int(cd8t.dict["migration"]):
                    shift = self.compute_shift(cd8t, distance, (nearest_tumour.xCOM, nearest_tumour.yCOM, nearest_tumour.zCOM))
                    self.move_cell(cd8t, shift)
                else:
                    
                    dx = nearest_tumour.xCOM - cd8t.xCOM
                    dy = nearest_tumour.yCOM - cd8t.yCOM
                    dz = nearest_tumour.zCOM - cd8t.zCOM
                    
                    end_x = cd8t.xCOM + int((dx/distance) * int(cd8t.dict["migration"]))
                    end_y = cd8t.yCOM + int((dy/distance) * int(cd8t.dict["migration"]))
                    end_z = cd8t.zCOM + int((dz/distance) * int(cd8t.dict["migration"]))
                    
                    shift = self.compute_shift(cd8t, int(cd8t.dict["migration"]), (end_x, end_y, end_z))
                    
                    self.move_cell(cd8t, shift)

               
class TumourCellsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        '''
        Establish a random direction for the tumour cell to move in.
        '''

        for tumour in self.cell_list_by_type(self.TUMOUR):
            tumour.lambdaVecX = tumour_migration * uniform(-0.5,0.5)
            tumour.lambdaVecY = tumour_migration * uniform(-0.5,0.5)


    def step(self, mcs):
        '''
        Move the tumour cell in a random direction at the tumour migration rate.
        '''
        
        for tumour in self.cell_list_by_type(self.TUMOUR):
            tumour.lambdaVecX = uniform(-0.5,0.5) * tumour_migration
            tumour.lambdaVecY = uniform(-0.5,0.5) * tumour_migration
            
class CAFsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        '''
        Establish a random direction for the CAF to move in.
        '''

        for caf in list(self.cell_list_by_type(self.CAF)) + list(self.cell_list_by_type(self.MYCAF)):
            caf.lambdaVecX = caf_migration * uniform(-0.5,0.5)
            caf.lambdaVecY = caf_migration * uniform(-0.5,0.5)


    def step(self, mcs):
        '''
        Move the CAF in a random direction at the CAF migration rate.
        '''
        
        for caf in list(self.cell_list_by_type(self.CAF)) + list(self.cell_list_by_type(self.MYCAF)):
            caf.lambdaVecX = caf_migration * uniform(-0.5,0.5)
            caf.lambdaVecY = caf_migration * uniform(-0.5,0.5)


########################################     
## CLASSES FOR CALIBRATING PARAMETERS ##
########################################
    
class CellMigrationTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
        self.step_counter = 0
        self.tumour_migrations = []
        self.caf_migrations = []
        
        self.cd8t_file_path = None
        self.cd8t_file_path = None
        self.cd8t_file_path = None
    #'''    
    def start(self):
        
        # Set up CSV file
        output_dir = self.output_dir
        self.cd8t_file_path = os.path.join(output_dir, "cd8t_migration.csv")
        self.tumour_file_path = os.path.join(output_dir, "tumour_migration.csv")
        self.caf_file_path = os.path.join(output_dir, "caf_migration.csv")
        
        with open(self.cd8t_file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MCS", "force", "migration"])
        f.close()
            
        with open(self.tumour_file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MCS", "migration"])
        f.close()
        
        with open(self.caf_file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MCS", "migration"])
        f.close()
        
        # Plot cd8t migration
        self.plot_cd8t_migration = self.add_new_plot_window(title="CD8 T migration",
                                                     x_axis_title='MonteCarlo Step (MCS)',
                                                     y_axis_title='migration per 10 MCS',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False,
                                                     config_options={'legend':True})
                                                     
        self.plot_cd8t_migration.add_plot("CD8 T migration",
            style='Lines', color='lightsteelblue', size=5) 
        
    def step(self, mcs):
        
        cd8t_migrations = []
        
        #default_cd8t_migration += 100
        #print(default_cd8t_migration)
                
        for cell in self.cell_list:
                
            dx = cell.xCOM - cell.dict["position_history"][0]
            dy = cell.yCOM - cell.dict["position_history"][1]
            dz = cell.zCOM - cell.dict["position_history"][2]
                
            displacement = math.sqrt(dx**2 + dy**2 + dz**2)
                                               
            cell.dict["position_history"][0] = cell.xCOM
            cell.dict["position_history"][1] = cell.yCOM
            cell.dict["position_history"][2] = cell.zCOM
                
            if cell.type == self.TUMOUR:
                if displacement < 5: # Remove artifact outliers
                    with open(self.tumour_file_path, "a", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow([mcs, (displacement)])
                    f.close()
                    self.tumour_migrations.append(displacement)                        
            elif cell.type == self.CAF or cell.type == self.MYCAF:
                with open(self.caf_file_path, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([mcs, (displacement)])
                f.close()
                self.caf_migrations.append(displacement)
            elif cell.type == self.CD8T:
                with open(self.cd8t_file_path, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([mcs, cell.dict["migration"], (displacement)])
                f.close()
                cd8t_migrations.append(displacement)
                    
                                    
        # Plot
                
        if len(cd8t_migrations) != 0:
            self.plot_cd8t_migration.add_data_point("CD8 T migration", mcs, sum(cd8t_migrations)/(len(cd8t_migrations)*10))
        else:
            self.plot_cd8t_migration.add_data_point("CD8 T migration", mcs, 0)
        
        self.step_counter += 1
    #'''
    
class CD8TKillAttemptsTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
        self.file_path = None
        
    def start(self):
        '''
        Set up CSV file for tracking CD8 T kill attempts.
        '''
        
        output_dir = self.output_dir
        self.file_path = os.path.join(output_dir, "cd8t_kill_attempts.csv")
        
        with open(self.file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Cell_ID", "Kill_Attempts", "End_Time (MCS)", "Kill_Attempts_/_MCS"])
        
    def finish(self):
        '''
        Calculate average rate of kill attempts by CD8 T cells.
        
        Returns:
            CSV file with the total number of kill attempts, lifespan, and kill rate per MCS of each CD8 T cell.
        '''
        
        # Store CD8 T kill attempt counts in CSV file
        
        with open(self.file_path, "a", newline="") as f:
            writer = csv.writer(f)
            for cell_id, info in self.shared_steppable_vars["cd8t_kill_attempts"].items():
                
                if info[1] == 0:
                    writer.writerow([cell_id, info[0], info[1], 0])
                else:
                    writer.writerow([cell_id, info[0], info[1], info[0]/info[1]])
        f.close()
        
    
#################################
## CLASSES FOR OUTPUTTING DATA ##
#################################

class OutputCSVSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.file_path = None
        
    def start(self): 
        '''
        Compute model outputs.
        
        Returns:
            CSV file for live and dead counts of each cell type over time.
            CSV file for count of cells expressing CD274, for each cell type.
            ...
        '''
        
        output_dir = self.output_dir
        self.cell_count_file_path = os.path.join(output_dir, "cell_count.csv")
        
        with open(self.cell_count_file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MCS", "Live Tumour", "Live CAF", "Live CD8 T", "Dead Tumour", "Dead CAF", "Dead CD8 T"])
                  
    def step(self, mcs):
        
        # Live cell count
        tumour_count = sum(1 for tumour in self.cell_list_by_type(self.TUMOUR))
        caf_count = sum(1 for caf in self.cell_list_by_type(self.CAF)) + sum(1 for mycaf in self.cell_list_by_type(self.MYCAF))
        cd8t_count = sum(1 for cd8t in self.cell_list_by_type(self.CD8T))
        
        with open(self.cell_count_file_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([mcs, tumour_count, caf_count, cd8t_count,
                self.shared_steppable_vars["dead_tumour_count"], self.shared_steppable_vars["dead_caf_count"],
                self.shared_steppable_vars["dead_cd8t_count"]])
        f.close()
        
        
            
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
        self.plot_live_count.add_plot("Live myCAF Count", style='Lines', color='gold', size=5)
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
  
        tumour_count = 0; caf_count = 0; mycaf_count = 0; cd8t_count = 0;
        
        # Plot live cell count
        for cell in self.cell_list:
                if cell.type == self. TUMOUR:
                    tumour_count += 1
                elif cell.type == self.CAF:
                    caf_count += 1
                elif cell.type == self.MYCAF:
                    mycaf_count += 1
                elif cell.type == self.CD8T:
                    cd8t_count += 1

        self.plot_live_count.add_data_point("Live Tumour Count", mcs, tumour_count)
        self.plot_live_count.add_data_point("Live CAF Count", mcs, caf_count)
        self.plot_live_count.add_data_point("Live myCAF Count", mcs, mycaf_count)
        self.plot_live_count.add_data_point("Live CD8T Count", mcs, cd8t_count)
        
        # Plot dead cell count
        self.plot_dead_count.add_data_point("Dead Tumour Count",
            mcs, self.shared_steppable_vars["dead_tumour_count"])
        self.plot_dead_count.add_data_point("Dead CD8T Count",
            mcs, self.shared_steppable_vars["dead_cd8t_count"])
        self.plot_dead_count.add_data_point("Dead CAF Count",
            mcs, self.shared_steppable_vars["dead_caf_count"])
            
          
        # Plot total diffusing agent
        
        # Fix to account for decay
        
        self.plot_diffusing_agent.add_data_point("Total IFN-Gamma", mcs, self.shared_steppable_vars["total_ifn_gamma"])
        self.plot_diffusing_agent.add_data_point("Total TGF-Beta", mcs, self.shared_steppable_vars["total_tgf_beta"])
        
        # Plot gene expression
        tumour_expressing_cd274 = sum(1 for tumour in self.cell_list_by_type(self.TUMOUR)
            if tumour.dict["CD274?"] == True)
        all_caf_expressing_cd274 = sum(1 for caf in self.cell_list_by_type(self.CAF) if caf.dict["CD274?"] == True) + \
            sum(1 for mycaf in self.cell_list_by_type(self.MYCAF) if mycaf.dict["CD274?"] == True)
        cd8t_expressing_cd274 = sum(1 for cd8t in self.cell_list_by_type(self.CD8T) 
            if cd8t.dict["CD274?"] == True)
        
        self.plot_gene_expression.add_data_point("Tumour With CD274",
            mcs, tumour_expressing_cd274)
        self.plot_gene_expression.add_data_point("CAF With CD274",
            mcs, all_caf_expressing_cd274)
        self.plot_gene_expression.add_data_point("CD8T With CD274",
            mcs, cd8t_expressing_cd274)
                                
        
        
      