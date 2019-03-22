'''
Created on Aug 15, 2016

@author: kfbz9246
'''

from itertools import product
from copy import deepcopy
from math import lgamma, exp

def interaction_test(data):
    marginals_table = {}
    order = len(data.columns)-1
    
    for face in product(['1', ':'], repeat=order):
        data_copy = data
        for column_index, column in enumerate(face):
            if column == '1':
                data_copy = data_copy[data_copy.iloc[:,column_index] == 0]
        marginals_table[face] = data_copy['Count'].sum()
    
    
    face_values_table = {}
    face_upper_left_is_added_table = {}
    upper_bound = 1000000000000000
    lower_bound = 0
    for face in product(['1', ':'], repeat=order):
        face_constraint_value = 0
        
        summed_column_count = 0
        for column in face:
            if column == ':':
                summed_column_count += 1
        for sub_face_component in product(['1', ':'], repeat=summed_column_count):
            sub_face_component = list(sub_face_component)
            summed_columns_replaced_count = 0
            next_sub_face = []
            for column in face:
                if column == ':':
                    next_sub_face_column = sub_face_component.pop(0)
                    if next_sub_face_column == '1':
                        summed_columns_replaced_count += 1
                    next_sub_face.append(next_sub_face_column)
                else:
                    next_sub_face.append(column)
            if summed_columns_replaced_count < summed_column_count:
                if summed_columns_replaced_count % 2 == 0:
                    face_constraint_value += marginals_table[tuple(next_sub_face)]
                else:
                    face_constraint_value -= marginals_table[tuple(next_sub_face)]
         
        if summed_column_count%2==0:
            face_values_table[face] = face_constraint_value
            face_upper_left_is_added_table[face] = True
            
            potential_lower_bound = -1*face_constraint_value
            if potential_lower_bound > lower_bound:
                lower_bound = potential_lower_bound
        else:
            face_values_table[face] = face_constraint_value
            face_upper_left_is_added_table[face] = False
            
            if face_constraint_value < upper_bound:
                upper_bound = face_constraint_value 
    print lower_bound
    print upper_bound
    
    test_range_lower_bound = lower_bound
    test_range_upper_bound = marginals_table[tuple(['1']*order)]

    if upper_bound-lower_bound >0:
        if float(marginals_table[tuple(['1']*order)]-lower_bound)/(upper_bound-lower_bound) >.5:
            test_range_lower_bound = marginals_table[tuple(['1']*order)]+1
            test_range_upper_bound = upper_bound

    p_value = 0
    for test_value in range(test_range_lower_bound, test_range_upper_bound+1):
        log_numerator_count = 0
        for cell_reference_part in product(['1', ':'], repeat=order-1):
            cell_reference_part = list(cell_reference_part)
            
            num_chosen = None
            num_chosen_cell_reference = tuple(['1']+cell_reference_part)
            num_chosen_from_cell_reference = tuple([':']+cell_reference_part)
            

            if face_upper_left_is_added_table[num_chosen_cell_reference]:
                num_chosen = face_values_table[num_chosen_cell_reference]+test_value
            else:
                num_chosen = face_values_table[num_chosen_cell_reference]-test_value
            
            num_chosen_from = None
            
    
            if face_upper_left_is_added_table[num_chosen_from_cell_reference]:
                num_chosen_from = face_values_table[num_chosen_from_cell_reference] + face_values_table[num_chosen_cell_reference] 
            else:
                num_chosen_from = face_values_table[num_chosen_from_cell_reference] + face_values_table[num_chosen_cell_reference]
                
            
            if num_chosen_from > 0:
                log_numerator_count += lgamma(num_chosen_from+1)-lgamma(num_chosen+1)-lgamma(num_chosen_from-num_chosen+1)
        numerator_count = exp(log_numerator_count)
            
        denominator_sum = 0.0
        for upper_left_corner_value in range(lower_bound, upper_bound+1):
            
            log_denominator_count = 0
            for cell_reference_part in product(['1', ':'], repeat=order-1):
                cell_reference_part = list(cell_reference_part)
                num_chosen_from = None
                num_chosen_cell_reference = tuple(['1']+cell_reference_part)
                num_chosen_from_cell_reference = tuple([':']+cell_reference_part)
                if face_upper_left_is_added_table[num_chosen_from_cell_reference]:
                    num_chosen_from = face_values_table[num_chosen_from_cell_reference] + face_values_table[num_chosen_cell_reference] 
                else:
                    num_chosen_from = face_values_table[num_chosen_from_cell_reference] + face_values_table[num_chosen_cell_reference]
                
                num_chosen = None
                
                if face_upper_left_is_added_table[num_chosen_cell_reference]:
                    num_chosen = face_values_table[num_chosen_cell_reference]+upper_left_corner_value
                else:
                    num_chosen = face_values_table[num_chosen_cell_reference]-upper_left_corner_value
                if num_chosen_from > 0:
                    log_denominator_count += lgamma(num_chosen_from+1)-lgamma(num_chosen+1)-lgamma(num_chosen_from-num_chosen+1)
    
            denominator_sum += exp(log_denominator_count)
        p_value += numerator_count/denominator_sum

    return p_value, upper_bound - lower_bound
    '''
    test_face = ('1','1','1')
    
    constraint_value = face_values_table[test_face]
    if face_upper_left_is_added_table[test_face]:
        print constraint_value + marginals_table[('1','1','1')]
    else:
        print constraint_value - marginals_table[('1','1','1')]
    '''
    

    
if __name__ == '__main__':
    from pandas import Series, DataFrame
    data = {'1' : Series([0,0,1,1]),
            '2' : Series([0,1,0,1]),
            'Count' : Series([100,100,100,100])}
    data = DataFrame(data)
    print data
    print interaction_test(data)


        
    