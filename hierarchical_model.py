'''
Created on Aug 25, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series
from numpy import mean, std, array, mean, var
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx
from itertools import combinations, izip
from test import interaction_test



microbe_genus = 'EC'
year = 2008
generate_new_graph_positions = False
show_order = False

engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year==year, Microbes.genus==microbe_genus)
suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year==year, Microbes.genus==microbe_genus)
data_query = resistant_query.union(suscptable_query)
data = read_sql(data_query.statement, engine)




#data.insert(2,'Count', 1)



data =  data.pivot('IsolateID', 'DrugID', 'Resistance')

data.dropna(1,thresh=int(len(data)*.2), inplace=True)
data.dropna(how='any',inplace=True)




data = data[['fox', 'tio', 'axo', 'amp', 'amc']]
print data

antibiotic_set = tuple(data.columns.values)

data.insert(0,'Count', 1)
data = data.groupby(antibiotic_set).sum().reset_index()


#graph = Graph()

for antibiotic_4_set in combinations(antibiotic_set,4):
    test_data = data[list(antibiotic_4_set)+['Count']]
    test_data = test_data.groupby(antibiotic_4_set).sum().reset_index()
    try:
        if interaction_test(test_data) <.05:
            all_three_way_interactions_signficant = True
            for antibiotic_triple in combinations(antibiotic_4_set,3):
                test_data_2 = data[list(antibiotic_triple)+['Count']]
                test_data_2 = test_data_2.groupby(antibiotic_triple).sum().reset_index()
                try:
                    if interaction_test(test_data_2) >.05:
                        all_three_way_interactions_signficant = False
                except OverflowError:
                    print 'three set untestable'
                    #all_three_way_interactions_signficant = False
            if all_three_way_interactions_signficant:
                print antibiotic_4_set
            #print test_data
            
    except OverflowError:
        print '4 set untestable'