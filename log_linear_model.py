'''
Created on Aug 24, 2016

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
from numpy.random import dirichlet
import statsmodels.formula.api as smf
from statsmodels.api import families

dimention = 3
antibiotic_set = ['chl', 'tet', 'cip']#['tet', 'fox', 'tio', 'axo', 'amp', 'amc']

microbe_genus = 'EC'

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


data =  data.pivot('IsolateID', 'DrugID', 'Resistance')

data.dropna(1,thresh=int(len(data)*.2), inplace=True)
data.dropna(how='any',inplace=True)





data = data[antibiotic_set]
antibiotic_set = tuple(data.columns.values)
data.insert(0,'Count', 1)
data = data.groupby(antibiotic_set).sum().reset_index()


rows = []
for row in product([0,1], repeat=len(antibiotic_set)):
    row_dict = {}
    for label, value in zip(antibiotic_set, row):
        row_dict[label] = value
    rows.append(row_dict)
    
all_resistance_combinations = DataFrame(rows)

data = merge(all_resistance_combinations, data, 'left', antibiotic_set)
data.fillna(0, inplace=True)

for antibiotic_pair in combinations(antibiotic_set, dimention):
    test_data = data[list(antibiotic_pair)+['Count']]
    test_data = test_data.groupby(antibiotic_pair).sum().reset_index()
    
    
    
    
    data_as_array = zeros([2]*dimention,dtype=int)
    for index, row in test_data.iterrows():
        data_as_array[tuple(row.iloc[range(dimention)])] = row['Count']
    print data_as_array
   
    all_columns_string = ''
    for column_name in antibiotic_pair:
        all_columns_string += column_name +'+'
    all_columns_string = all_columns_string[0:-1]
    
    model = smf.glm(formula="Count ~ ("+ all_columns_string +")**3-1", data=test_data, family=families.Poisson())
    model_fit = model.fit()
    print model_fit.summary()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    