'''
Created on Sep 15, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, DataFrame, merge
from numpy import mean, std, array, mean, var, zeros
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log10, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx
from itertools import combinations, product
from test import interaction_test
from _sqlite3 import Row
from numpy.random import dirichlet



antibiotic_set = ['axo', 'tio', 'amp']#['chl', 'tet', 'cip']#

microbe_genus = 'EC'

generate_new_graph_positions = False
show_order = False

engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=2011, Samples.year<=2013, Microbes.genus==microbe_genus)
suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=2011, Samples.year<=2013, Microbes.genus==microbe_genus)
data_query = resistant_query.union(suscptable_query)
data = read_sql(data_query.statement, engine)


data =  data.pivot('IsolateID', 'DrugID', 'Resistance')

data.dropna(1,thresh=int(len(data)*.2), inplace=True)
data.dropna(how='any',inplace=True)




'''
rows = []
for row in product([0,1], repeat=len(antibiotic_set)):
    row_dict = {}
    for label, value in zip(antibiotic_set, row):
        row_dict[label] = value
    rows.append(row_dict)
    
all_resistance_combinations = DataFrame(rows)

data = merge(all_resistance_combinations, data, 'left', antibiotic_set)
data.fillna(0, inplace=True)
'''
#data.sort_values([item for item in reversed(antibiotic_set)], inplace=True)
for antibiotic_subset in combinations(antibiotic_set, 3):
    antibiotic_subset = list(antibiotic_subset)
    data_subset = data[antibiotic_subset]
    antibiotic_set = tuple(data_subset.columns.values)
    data_subset.insert(0,'Count', 1)
    data_subset = data_subset.groupby(antibiotic_set).sum().reset_index()
    data_subset.sort_values(antibiotic_subset, inplace=True)
    
    data_as_array = zeros([2]*3,dtype=int)
    for index, row in data_subset.iterrows():
        data_as_array[tuple(row.iloc[range(3)])] = row['Count']
    print antibiotic_subset
    print data_as_array
    
    print data_subset[antibiotic_subset+ ['Count']]

#data[antibiotic_set+ ['Count']].to_csv('drug_table.csv')