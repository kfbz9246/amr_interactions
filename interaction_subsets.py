'''
Created on Sep 8, 2017

@author: root
'''
from pandas import read_excel, Series, read_csv

four_way_interactoins = read_csv('4.csv', header=0, index_col=0)
three_way_interactions = read_csv('3.csv', header=0, index_col=0)


for significance_column in ['p_value_Chickens', 'p_value_Turkeys']:
    four_way_sets = []
    for name, data in four_way_interactoins[significance_column].iteritems():
        if data < .01:
            four_way_sets.append(set(name.split(':')))
    
    three_way_set_to_name_dict = {}
    three_way_sets = []
    for name, data in three_way_interactions[significance_column].iteritems():
        if data < .01:
            new_three_way_set = set(name.split(':'))
            three_way_sets.append(new_three_way_set)
            three_way_set_to_name_dict[str(new_three_way_set)] = name
    
    covered_3_sets = set()
    for four_way_set in four_way_sets:
        for three_way_set in three_way_sets:
            if three_way_set.issubset(four_way_set):
                covered_3_sets.add(three_way_set_to_name_dict[str(three_way_set)])
    
    three_way_set_is_covered_dict = {}
    for covered_3_set in covered_3_sets:
        three_way_set_is_covered_dict[covered_3_set] = 1
    three_way_interactions[significance_column+'_covered'] = Series(three_way_set_is_covered_dict)
three_way_interactions.to_csv('3_wauy_interactions_coverage.csv')
