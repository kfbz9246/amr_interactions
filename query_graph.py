'''
Created on Oct 20, 2016

@author: kfbz9246
'''
from pickle import load
from itertools import combinations

for line in open('foo.txt'):
    print line
    antibiotic_set = line.split()
    
    
    graph = load(open('graph.pickle'))
    is_clique = True
    for combination in combinations(antibiotic_set, 2):
        if not graph.has_edge(*combination):
            is_clique =False
            
    print str(antibiotic_set) + '-' + str(is_clique)