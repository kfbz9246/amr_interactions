'''
Created on Aug 24, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series
from numpy import mean, std, array, mean, var, zeros
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx
from itertools import combinations, izip
from test import interaction_test


dimention = 4
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
antibiotic_set = tuple(data.columns.values)

data.insert(0,'Count', 1)
data = data.groupby(antibiotic_set).sum().reset_index()


number_tables_list = {}
for antibiotic_pair in combinations(antibiotic_set, dimention):
    test_data = data[list(antibiotic_pair)+['Count']]
    test_data = test_data.groupby(antibiotic_pair).sum().reset_index()
    data_as_array = zeros([2]*dimention,dtype=int)
    for index, row in test_data.iterrows():
        data_as_array[tuple(row.iloc[range(dimention)])] = row['Count']
    print data_as_array
    try:
        p_value, number_viable_tables = interaction_test(test_data) 
        if not number_tables_list.has_key(number_viable_tables):
            number_tables_list[number_viable_tables] = 0
        number_tables_list[number_viable_tables] += 1
        print p_value
        
            
    except OverflowError:
        print 'overflow'

temp_fig = figure(figsize=(9.25, 7), dpi=96, frameon=True)
fig_ax = temp_fig.add_subplot(111)

min_table_size = 10000
max_table_size = 1
for number_viable_tables, count in number_tables_list.iteritems():
    if number_viable_tables + 1 > max_table_size:
        max_table_size = number_viable_tables + 1
    if number_viable_tables +1 < min_table_size:
        min_table_size = number_viable_tables + 1
    for index in range(count):
        fig_ax.scatter([number_viable_tables+1], [index+1])
fig_ax.set_title('Number of Tables with Given Margins\n'+ str(dimention) +'-Way Tables')
fig_ax.set_xlabel('Number of Tables')
fig_ax.set_ylabel('Count')
fig_ax.set_xlim(min_table_size-1, max_table_size+1)
savefig('figure_1.png', bbox_inches='tight')
