'''
Created on Jul 18, 2016

@author: kfbz9246
'''
from scipy.misc import comb
    
class ContingencyValue(object):
    def __init__(self, is_upper_left_corner, added, indices):
        self.is_upper_left_corner = is_upper_left_corner
        self.added = added
        self.indicies = indices
        self.value = None
        
    def is_added(self):
        return self.added
    
        
    def split(self):
        child_1 = [ContingencyValue(self.is_upper_left_corner, self.added, self.indicies+[0])]
        child_2 = None
        if self.is_upper_left_corner:
            child_2_1 = ContingencyValue(False, self.added, self.indicies+[':'])
            if self.added:
                new_added = False
            else:
                new_added = True   
            child_2_2 = ContingencyValue(True, new_added, self.indicies+['0'])
            child_2 = [child_2_1, child_2_2]
        else:
            child_2 =[ContingencyValue(self.is_upper_left_corner, self.added, self.indicies+[1])]  
        return child_1, child_2
    
    
    def is_upper_left_corner(self):
        return self.is_upper_left_corner
    
    def instantiate(self, marginals_table):
        if not self.is_upper_left_corner:
            sub_contigency_table = marginals_table
            for data_index, contigency_value_index in enumerate(self.indices):
                if contigency_value_index != ':':
                    sub_contigency_table = sub_contigency_table[sub_contigency_table.iloc[data_index] == contigency_value_index]
            self_value = sum(sub_contigency_table['count'])
            if not self.added:
                self_value =-1*self_value
            self.value = self_value

    def evaluate(self, upper_left_corner):
        return_val = upper_left_corner
        if not self.is_upper_left_corner:
            return_val = self.value
        return return_val
        
    def __str__(self):
        return_string = 'n'
        for contingency_table_index in self.indicies:
            return_string += str(contingency_table_index)
        return return_string
            
class Combination(object):
    def __init__(self, choosing_from=ContingencyValue(False, True, [':']), chosen=[ContingencyValue(True,True, [0])]):
        self.choosing_from = choosing_from
        self.chosen = chosen
    
    def split(self):
        new_choosing_from = self.choosing_from.split()
        new_choosen_1 = []
        new_choosen_2 = []
        for chosen_value in self.chosen:
            chosen_value_split = chosen_value.split()
            new_choosen_1.extend(chosen_value_split[0])
            new_choosen_2.extend(chosen_value_split[1])
        child_1 = Combination(new_choosing_from[0][0], new_choosen_1)
        child_2 = Combination(new_choosing_from[1][0], new_choosen_2)
        return child_1, child_2
    
    def instantiate(self, marginals_table):
        self.choosing_from.instantiate(marginals_table)
        for contingency_value in self.chosen:
            contingency_value.instantiate(marginals_table)

    def evaluate(self, upper_left_corner_value):
        choosing_from_val = self.choosing_from.evaluate(upper_left_corner_value)
        total_chosen = 0
        for chosen_value in self.chosen:
            total_chosen += chosen_value.evaluate(upper_left_corner_value)
        return comb(choosing_from_val, total_chosen)

    def __str__(self):
        return_string = 'comb(' + str(self.choosing_from) + ', '
        for chosen_value in self.chosen:
            if chosen_value.is_added():
                return_string += '+'
            else:
                return_string += '-'
            return_string += str(chosen_value)
        return_string += ')'
        return return_string
    
class Formula(object):
    def __init__(self):
        self.combination_list = [Combination()]
    
    def split(self):
        combination_list_1 = []
        combination_list_2 = []
        for combination in self.combination_list:
            combination_split = combination.split()
            combination_list_1.append(combination_split[0])
            combination_list_2.append(combination_split[1])
        self.combination_list = combination_list_1 + combination_list_2
    
    def instantiate(self, marginals_table):
        for combination in self.combination_list:
            combination.instantiate(marginals_table)
    
    def evaluate(self, upper_left_corner_value):
        return_val = 1
        for combination in self.combination_list:
            return_val *= combination.evaluate(upper_left_corner_value)
        return return_val
    
    def __str__(self):
        return_string = ''
        for combination in self.combination_list:
            return_string += str(combination)
            return_string += '*'
        return return_string

if __name__ == '__main__':
    my_formula = Formula()
    split_val = 1
    for i in range(split_val):
        my_formula.split()
    print str(my_formula)
    
