'''
Created on Jul 26, 2016

@author: kfbz9246
'''
from itertools import product, izip_longest
from copy import deepcopy

order = 3

class Constraint(object):
    def __init__(self, cell_reference, margin_index):
        self.cell_reference = cell_reference
        self.margin = deepcopy(cell_reference)
        self.margin[margin_index] = ':'
        left_hand_side = ConstraintChunk(True, cross_product).resolve()
        for item in left_hand_side:
            item.switch_is_added()
        right_hand_side = ConstraintChunk(True, self.margin).resolve()
        self.constraint = right_hand_side + left_hand_side
        self.constraint.sort()
        constraint_index = 0
        while constraint_index < len(self.constraint)-1:
            if self.constraint[constraint_index] == self.constraint[constraint_index+1] and self.constraint[constraint_index].get_is_added() != self.constraint[constraint_index+1].get_is_added():
                    self.constraint.pop(constraint_index + 1)
                    self.constraint.pop(constraint_index)
            else:
                constraint_index += 1
        
        
    def get_cell_reference(self):
        return self.cell_reference
    
    def get_margin(self):
        return self.margin
    
    def get_constraint(self):
        return self.constraint
    
    def __cmp__(self, other):
        cmp_val = 0
        for constraint_chunk_1, constraint_chunk_2 in izip_longest(self.constraint, other.constraint):
            cmp_val =  cmp(constraint_chunk_1, constraint_chunk_2)
            if cmp_val != 0:
                break
        return cmp_val
    
class ConstraintChunk(object):
    def __init__(self, is_added, cell_index):
        self.is_added = is_added
        self.cell_index = cell_index
        self.upper_left_corner = True
        for dimention in self.cell_index:
            if dimention != 1:
                self.upper_left_corner = False
         
    def resolve(self):
        resolved_indices_left = []
        resolved_indices_right = []
        two_not_found = True
        for dimention in self.cell_index:
            if dimention == 2 and two_not_found:
                two_not_found = False
                resolved_indices_left.append(':')
                resolved_indices_right.append(1)
            else:
                resolved_indices_left.append(dimention)
                resolved_indices_right.append(dimention)
        return_val = []
        if two_not_found:
            return_val.append(self)
        else:
            return_val.extend(ConstraintChunk(self.is_added,resolved_indices_left).resolve())
            return_val.extend(ConstraintChunk(False if self.is_added else True, resolved_indices_right).resolve())
        return return_val
    
    def switch_is_added(self):
        if self.is_added:
            self.is_added = False
        else:
            self.is_added = True
            
    def get_is_added(self):
        return self.is_added
    
    def is_upper_left_corner(self):
        return self.upper_left_corner
                
    def __cmp__(self, other):
        cmp_val = None
        if type(self) != type(other):
            cmp_val = -1
        else:
            for  self_dimention, other_dimention in zip(self.cell_index, other.cell_index):
                cmp_val = cmp(str(self_dimention), str(other_dimention))
                if cmp_val != 0:
                    break
        return cmp_val
    
    def __str__(self):
        return_val = None
        if self.is_added:
            return_val='{'+str(self.cell_index)+'}'
        else:
            return_val = '{-'+str(self.cell_index)+'}'
        return return_val
        
    def __repr__(self):
        return str(self)
    
constraint_list = []
for cross_product in product([1,2], repeat=order):
    cross_product = list(cross_product)
    for cross_product_index in range(len(cross_product)):
        constraint_list.append(Constraint(cross_product, cross_product_index))

        

constraint_list_index = 0
while constraint_list_index < len(constraint_list):
    next_constraint = constraint_list[constraint_list_index]
    equivelent_constraints = [constraint_list[constraint_list_index]]
    constraint_list_index_2 = constraint_list_index + 1
    while constraint_list_index_2 < len(constraint_list):
        if cmp(constraint_list[constraint_list_index], constraint_list[constraint_list_index_2]) == 0:
            equivelent_constraints.append(constraint_list.pop(constraint_list_index_2))
        else:
            constraint_list_index_2 += 1
    
    #print equivelent_constraints[0].get_constraint()
    for constraint in equivelent_constraints:
        print str(constraint.get_cell_reference()) + ' | ' + str(constraint.get_margin()) + ' | ' + str(constraint.get_constraint()) + ' | ' + str(len(constraint.get_constraint()))
    
    print'---'
    constraint_list_index += 1

'''for constraint in constraint_list:
    print constraint

print len(constraint_list)
'''         
            
    