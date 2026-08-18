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

# Set established parameters (will not change)

tumour_lambda_vol = 50 
caf_lambda_vol = 100
cd8t_lambda_vol = 100


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
                    cell.targetVolume = random.uniform(0.25, 1.75) * caf_vol
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
            CSV file with all outputs per MCS.
        '''
     
        # Initialize lists of variables
                
        self.shared_steppable_vars["cd8t_kill_attempts_set"] = {}
        total_mcs = self.simulator.getNumSteps()
        for cd8t in self.cell_list_by_type(self.CD8T):
            self.shared_steppable_vars["cd8t_kill_attempts_set"][cd8t.id] = [0, total_mcs]
            
        self.shared_steppable_vars["caf_density_list"] = []            
        self.shared_steppable_vars["cd8t_tumour_dist_list"] = []
        self.shared_steppable_vars["cd8t_caf_dist_list"] = []
        self.shared_steppable_vars["cd8t_density_list"] = []
                        
        # Initialize shared variables
            
        self.shared_steppable_vars["initial_tumour_count"] = sum(1 for tumour in self.cell_list_by_type(self.TUMOUR)) 
        self.shared_steppable_vars["tumour_divisions"] = 0
        self.shared_steppable_vars["tumour_apoptosis"] = 0
        
        self.shared_steppable_vars["caf_count"] = 0
        self.shared_steppable_vars["caf_density"] = 0
        self.shared_steppable_vars["caf_divisions"] = 0
        
        self.shared_steppable_vars["cd8t_kill_attempts"] = 0
        self.shared_steppable_vars["cd8t_successful_kills"] = 0
        self.shared_steppable_vars["cd8t_apoptosis"] = 0
              
               
        # Set up CSV file
        
        output_dir = self.output_dir
        self.output_file_path = os.path.join(output_dir, "output.csv")
        
        with open(self.output_file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MCS",
                "Tumour Fold Change", "Tumour Total Divisions", "Tumour Total Natural Apoptosis",
                "CAF Count", "Local CAF Density", "CAF Total Divisions",
                "CD8 T Count", "Exhausted CD8 T Count", "CD8 T Kill Attempts", "CD8 T Successful Kills",
                "CD8 T Kill Efficiency", "CD8 T Apoptosis", "Mean Nearest CD8 T Tumour Distance", "Mean Nearest CD8 T CAF Distance",
                "Mean Local CD8 T Density Around Tumour"])
        f.close()
            
    def step(self, mcs):
        
        # Skip first MCS
        if mcs == 0:
            return
        
        # CALCULATE VARIABLES
        
        # Tumour fold change
        if self.shared_steppable_vars["initial_tumour_count"] == 0:
            tumour_fold_change = 1
        else:
            tumour_fold_change = (sum(1 for tumour in self.cell_list_by_type(self.TUMOUR))) / self.shared_steppable_vars["initial_tumour_count"]
        
        # CAF cell count
        caf_count = sum(1 for caf in self.cell_list_by_type(self.CAF)) + sum(1 for mycaf in self.cell_list_by_type(self.MYCAF))
        
        # Mean local CAF density around tumour
        if len(self.shared_steppable_vars["caf_density_list"]) == 0:
            caf_density = 0
        else:
            caf_density = sum(self.shared_steppable_vars["caf_density_list"]) / len(self.shared_steppable_vars["caf_density_list"])
                        
        # CD8 T cell count
        cd8t_count = sum(1 for cd8t in self.cell_list_by_type(self.CD8T))
        
        # Exhausted CD8 T cell count
        exhausted_cd8t_count = sum(1 for cd8t in self.cell_list_by_type(self.CD8T) if cd8t.dict["exhaustion_counter"] >= exhaustion_threshold)
        
        # CD8 T kill efficiency
        if self.shared_steppable_vars["cd8t_kill_attempts"] == 0:
            cd8t_kill_efficiency = 0
        else:
            cd8t_kill_efficiency = self.shared_steppable_vars["cd8t_successful_kills"] / self.shared_steppable_vars["cd8t_kill_attempts"]
        
        # Mean nearest CD8 T Tumour distance
        if len(self.shared_steppable_vars["cd8t_tumour_dist_list"]) == 0:
            cd8t_tumour_dist = 0
        else:
            cd8t_tumour_dist = sum(self.shared_steppable_vars["cd8t_tumour_dist_list"]) / len(self.shared_steppable_vars["cd8t_tumour_dist_list"])
        
        # Mean nearest CD8 T CAF distance
        if len(self.shared_steppable_vars["cd8t_caf_dist_list"]) == 0:
            cd8t_caf_dist = 0
        else:
            cd8t_caf_dist = sum(self.shared_steppable_vars["cd8t_caf_dist_list"]) / len(self.shared_steppable_vars["cd8t_caf_dist_list"])
        
        # Mean local CD8 T density around tumour
        if len(self.shared_steppable_vars["cd8t_density_list"]) == 0:
            cd8t_density = 0
        else:
            cd8t_density = sum(self.shared_steppable_vars["cd8t_density_list"]) / len(self.shared_steppable_vars["cd8t_density_list"])
        
        
        # Write to CSV file
        
        with open(self.output_file_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([mcs, 
                tumour_fold_change, self.shared_steppable_vars["tumour_divisions"], self.shared_steppable_vars["tumour_apoptosis"],
                caf_count, caf_density, self.shared_steppable_vars["caf_divisions"],
                cd8t_count, exhausted_cd8t_count, self.shared_steppable_vars["cd8t_kill_attempts"], self.shared_steppable_vars["cd8t_successful_kills"],
                cd8t_kill_efficiency, self.shared_steppable_vars["cd8t_apoptosis"], cd8t_tumour_dist, cd8t_caf_dist,
                cd8t_density])
        f.close()
        
        # Clear lists for next MCS
        self.shared_steppable_vars["cd8t_tumour_dist_list"] = []
        self.shared_steppable_vars["cd8t_caf_dist_list"] = []
        

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
                cell.targetVolume *= 1 + tumour_prolif
            elif cell.type == self.CAF or cell.type == self.MYCAF:
                cell.targetVolume *= 1 + caf_prolif    

        
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
                    self.shared_steppable_vars["tumour_divisions"] += 1
                    cells_to_divide.append(cell)
            elif cell.type == self.CAF or cell.type == self.MYCAF:
                if cell.volume >= 2*caf_vol:
                    self.shared_steppable_vars["caf_divisions"] += 1
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
            
        
    def step(self, mcs):
        '''
        Update attributes of tumour cells and their effects on other cells.
        '''
        
        cells_to_delete = []
        
        # Track outputs
        
        # Local CD8 T density
        cd8t_cells = list(self.cell_list_by_type(self.CD8T)) 
        
        if len(cd8t_cells) != 0:      
            cd8t_positions = [(cd8t.xCOM, cd8t.yCOM, cd8t.zCOM) for cd8t in cd8t_cells]
            cd8t_tree = KDTree(cd8t_positions)
            
        # Local CAF density
        cafs = list(self.cell_list_by_type(self.CAF)) + list(self.cell_list_by_type(self.MYCAF))
        
        if len(cafs) != 0:      
            caf_positions = [(caf.xCOM, caf.yCOM, caf.zCOM) for caf in cafs]
            caf_tree = KDTree(caf_positions)
        
        
        for i, tumour in enumerate(self.cell_list_by_type(self.TUMOUR)):
            
            # Track outputs
            
            # Local CD8 T density
            if len(cd8t_cells) == 0:
                self.shared_steppable_vars["cd8t_density_list"].append(0)
            else:
                indices = cd8t_tree.query_ball_point((tumour.xCOM, tumour.yCOM, tumour.zCOM), r=10)
                self.shared_steppable_vars["cd8t_density_list"].append(len(indices))
                
            # Local CAF density
            if len(cafs) == 0:
                self.shared_steppable_vars["caf_density_list"].append(0)
            else:
                indices = caf_tree.query_ball_point((tumour.xCOM, tumour.yCOM, tumour.zCOM), r=10)
                self.shared_steppable_vars["caf_density_list"].append(len(indices))
            
            # Baseline apoptosis
            if random.random() <= tumour_apoptosis_prob:
                self.shared_steppable_vars["tumour_apoptosis"] += 1
                cells_to_delete.append(tumour)
                continue

            # CHECK 1: Done in CD8TCellsMoveSteppable                                
                    
            # CHECK 2: if CD274 expression is induced
            if self.field.IFN_gamma[int(tumour.xCOM), int(tumour.yCOM), int(tumour.zCOM)] > tumour_ifn_pdl1_threshold:
                tumour.dict["CD274?"] = True
            
            self.helper_func.update_lattice_sites(tumour.xCOM, tumour.yCOM, tumour.zCOM, self.field.TGF_beta,
                tumour.targetVolume, tumour_tgf_secretion)
        
        # Delete tumour cells marked for apoptosis die
        for tumour in cells_to_delete:
            self.delete_cell(tumour)


class UpdateCAFsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
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

            x1 = -0.00000000002*(tgf**2) + 449882*tgf + 0.1084
            x2 = 1 - x1
            
            if x2 <= 0:
              prob = 1
            else:
              prob = 1 - x2 ** (1/2160)
            
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
        
        # Delete CAFs marked for apoptosis die
        for caf in cells_to_delete:
            self.delete_cell(caf)
        

class UpdateCD8TCellsSteppable(SteppableBasePy):
    
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.helper_func = HelperFunctionsSteppable()
               
    def step(self, mcs):
        '''
        Update attributes of CD8 T cells and their effects on other cells.
        '''
        
        cells_to_delete = []
        
        cafs = list(self.cell_list_by_type(self.CAF)) + list(self.cell_list_by_type(self.MYCAF))
        
        if len(cafs) != 0:      
            caf_positions = [(caf.xCOM, caf.yCOM, caf.zCOM) for caf in cafs]
            caf_tree = KDTree(caf_positions)
        
        for cd8t in self.cell_list_by_type(self.CD8T):
            
            # Compute mean nearest distance between CD8 T cells & CAFs
            
            if len(cafs) != 0:
                 
                distance, index = caf_tree.query((cd8t.xCOM, cd8t.yCOM, cd8t.zCOM))
                self.shared_steppable_vars["cd8t_caf_dist_list"].append(distance)
                       
            # Apoptosis rate
            if random.random() <= cd8t_apoptosis_prob:
                self.shared_steppable_vars["cd8t_apoptosis"] += 1
                cells_to_delete.append(cd8t)
                continue
                
            # Check 1 done in UpdateTumourCellsSteppable and UpdateCAFsSteppable
            # Check 2 done in UpdateTumourCellsSteppable
                
            # CD8 T cell migration is affected by collagen density           
            collagen = int(self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM])
            
            cd8t.dict["migration"] = default_cd8t_migration
            
                        
            collagen = self.field.Collagen[cd8t.xCOM, cd8t.yCOM, cd8t.zCOM]
            if collagen <= 6.66034:
                cd8t.dict["migration"] = -1.1376 * collagen + 7.5768
            else:
                cd8t.dict["migration"] = 0
            
        # Delete CD8 T cells marked for apoptosis 
        for cd8t in cells_to_delete:
            
            # Track CD8 T kills
            self.shared_steppable_vars["cd8t_kill_attempts_set"][cd8t.id][1] = mcs
            
            self.delete_cell(cd8t)
            
            
###############################
## CLASSES FOR CELL MOVEMENT ##
###############################


class CD8TCellsMoveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.helper_func = HelperFunctionsSteppable()
        
    def can_shift(self, cell, shift):
        '''
        Check if all the cell's pixels can move by the shift
        '''
        
        dims = self.cell_field.getDim()
        
        for pt in self.get_cell_pixel_list(cell):
            x = pt.pixel.x + shift[0]
            y = pt.pixel.y + shift[1]
            z = pt.pixel.z + shift[2]
            
            if not (0 <= x < dims.x and 0 <= y < dims.y and 0 <= z < dims.z):
                return False
                
            occupant = self.cell_field[x, y, z]
            if occupant is not None and occupant.id != cell.id:
                return False
        
        # extra buffer only against same-type cells
        '''
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                nx, ny = x + dx, y + dy
                if 0 <= nx < dims.x and 0 <= ny < dims.y:
                    neighbor = self.cell_field[nx, ny, z]
                    if neighbor is not None and neighbor.type == cell.type and neighbor.id != cell.id:
                        return False
        '''
        return True

        
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
        
        for i in range (1, int(distance) + 1):
            t = i / distance
            x = int(x0 + (x1 - x0) * t)
            y = int(y0 + (y1 - y0) * t)
            z = int(z0 + (z1 - z0) * t)
            
            shift = (int(round(x - x0)), int(round(y - y0)), int(round(z - z0)))
        
            if self.can_shift(cell, shift):
                dest_x = x; dest_y = y; dest_z = z
            else:
                break
                
        return (int(round(dest_x - x0)), int(round(dest_y - y0)), int(round(dest_z - z0)))


    def step(self, mcs):
        
                 
        tumour_cells = list(self.cell_list_by_type(self.TUMOUR))
        
        if len(tumour_cells) == 0:
            return
        
        tumour_positions = [(tumour.xCOM, tumour.yCOM, tumour.zCOM) for tumour in tumour_cells ]
        tumour_tree = KDTree(tumour_positions)
        
        tumours_to_kill = []
        cell_ids_to_delete = set()
         
        for cd8t in self.cell_list_by_type(self.CD8T):
            
            distance, index = tumour_tree.query((cd8t.xCOM, cd8t.yCOM, cd8t.zCOM))
            nearest_tumour = tumour_cells[index]
                        
            self.shared_steppable_vars["cd8t_tumour_dist_list"].append(distance) # Track outputs
                                
            if distance > 0:     
                                              
                # Shift
                if distance <= int(cd8t.dict["migration"]):
                    if nearest_tumour not in tumours_to_kill and nearest_tumour.id not in cell_ids_to_delete:
                        tumours_to_kill.append(nearest_tumour)
                        shift = self.compute_shift(cd8t, distance, (nearest_tumour.xCOM, nearest_tumour.yCOM, nearest_tumour.zCOM))
                        self.move_cell(cd8t, shift)
                    else:
                        print("Avoid CD8 T cell collision")
                else:
                    
                    dx = nearest_tumour.xCOM - cd8t.xCOM
                    dy = nearest_tumour.yCOM - cd8t.yCOM
                    dz = nearest_tumour.zCOM - cd8t.zCOM
                    
                    end_x = cd8t.xCOM + int((dx/distance) * int(cd8t.dict["migration"]))
                    end_y = cd8t.yCOM + int((dy/distance) * int(cd8t.dict["migration"]))
                    end_z = cd8t.zCOM + int((dz/distance) * int(cd8t.dict["migration"]))
                    
                    shift = self.compute_shift(cd8t, int(cd8t.dict["migration"]), (end_x, end_y, end_z))
                    
                    if abs(shift[0] - dx) <= 3 and abs(shift[1] - dy) <=3 and abs(shift[2] - dz) <=3:
                        if nearest_tumour not in tumours_to_kill  and nearest_tumour.id not in cell_ids_to_delete:
                            tumours_to_kill.append(nearest_tumour)
                            self.move_cell(cd8t, shift)
                        else:
                            print("Avoid CD8 T cell collision")
                    else:
                        self.move_cell(cd8t, shift)  
                        
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cd8t):
                if neighbor:
                    if neighbor.type == self.TUMOUR:

                        tumour = neighbor
                        
                        # CD8 T cell secretes IFN gamma
                        self.helper_func.update_lattice_sites(cd8t.xCOM, cd8t.yCOM, cd8t.zCOM, self.field.IFN_gamma,
                            tumour.targetVolume, cd8t_ifn_secretion)
                        
                        # Track CD8 T kills
                        self.shared_steppable_vars["cd8t_kill_attempts"] += 1
                        self.shared_steppable_vars["cd8t_kill_attempts_set"][cd8t.id][0] += 1
                        
                        # Check if CD8 T cell successfully kills tumour cell:
                        if (cd8t.dict["CD274?"] == False or tumour.dict["CD274?"] == False) and cd8t.dict["exhaustion_counter"] < exhaustion_threshold:
                            cd8t.dict["exhaustion_counter"] += 1
                            self.shared_steppable_vars["cd8t_successful_kills"] += 1
                            cell_ids_to_delete.add(tumour.id)
                            continue
                        # Else, immune escape and CD8 T exhaustion occurs
                        else:
                            cd8t.dict["exhaustion_counter"] = exhaustion_threshold
                  
        for cell_id in cell_ids_to_delete:
          cell = self.inventory.attemptFetchingCellById(cell_id)
          if cell is not None:
            self.delete_cell(cell)
          else:
            print("CELL NO LONGER EXISTS WOMP WOMP")
    #'''           
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
        
    # Used to determine cell migration parameters    
    '''    
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
            for cell_id, info in self.shared_steppable_vars["cd8t_kill_attempts_set"].items():
                
                if info[1] == 0:
                    writer.writerow([cell_id, info[0], info[1], 0])
                else:
                    writer.writerow([cell_id, info[0], info[1], info[0]/info[1]])
        f.close()
        
        
      
