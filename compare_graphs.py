'''
Created on Jul 30, 2017

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame
from numpy import std, array, mean, var, zeros
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log10, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx, write_graphml, adjacency_matrix
from itertools import combinations, izip, product
from numpy.random import dirichlet
import statsmodels.formula.api as smf
from statsmodels.api import families

antibiotic_query_set = []#['tio','amc','axo','fox']#['amc', 'fox', 'axo', 'tio']#['amc', 'axo', 'fox']#['amp', 'tio', 'amc', 'fox', 'axo', 'chl']#['cot', 'gen', 'str', 'fis', 'tet', 'kan', 'chl', 'amp'] #['cip', 'gen', 'fis']#['chl', 'tet', 'cip']#
dimention = 2
microbe_genus = 'EC'
host = 'Chickens'
start_year = 2011
end_year = 2013


number_simulations = 1000000
alpha = .00001




graphs_dict = {}
source_list = ['Slaughter', 'Retail', 'All']
for source in source_list:
    print source
    engine = create_engine('sqlite:////home/workspace/amr_db/narms_database.sqlite')
    session = sessionmaker(bind=engine)()
    
    resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                         join(Isolates, Microbes, Samples).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
    suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                         join(Isolates, Microbes, Samples).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
    
    if len(antibiotic_query_set) > 0:
        resistant_query = resistant_query.filter(Tests.drug_id.in_(antibiotic_query_set))
        suscptable_query = suscptable_query.filter(Tests.drug_id.in_(antibiotic_query_set))
    
    if source != 'All':
        resistant_query = resistant_query.filter(Samples.source == source)
        suscptable_query = suscptable_query.filter(Samples.source == source)
        
    
    data_query = resistant_query.union(suscptable_query)
    data = read_sql(data_query.statement, engine)
    
     
    
    data =  data.pivot('IsolateID', 'DrugID', 'Resistance')
    
    data.dropna(1,thresh=int(len(data)*.2), inplace=True)
    data.dropna(how='any',inplace=True)
    antibiotic_set = antibiotic_query_set
    if len(antibiotic_query_set) == 0:
        antibiotic_set = tuple(data.columns.values)
        print antibiotic_set
    
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
    
    
    graph = Graph()
    for antibiotic in antibiotic_set:
        graph.add_node(antibiotic)
    
    
    for antibiotic_pair in combinations(antibiotic_set, dimention):
        test_data = data[list(antibiotic_pair)+['Count']]
        test_data = test_data.groupby(antibiotic_pair).sum().reset_index()
    
        cell_counts_list = []
        cell_counts_even = []
        for row_index, row in test_data.iterrows():
            cell_index_sum = 0
            for element_index, element in row.iteritems():
                if element_index == 'Count':
                    cell_counts_list.append(element+1)
                else:
                    cell_index_sum += element
            if cell_index_sum % 2==0:
                cell_counts_even.append(True)
            else:
                cell_counts_even.append(False)
        
    
        raw_odds_ratio = 1
        for is_even, sample_value in zip(cell_counts_even, cell_counts_list):
            if is_even:
                raw_odds_ratio *= sample_value-.5
            else:
                raw_odds_ratio /= sample_value-.5
        raw_odds_ratio = log10(raw_odds_ratio)
        
        sample_odds_ratios = []
        for i in range(number_simulations):
            sample = dirichlet(cell_counts_list)
            odds_ratio = 1
            for is_even, sample_value in zip(cell_counts_even, sample):
                if is_even:
                    odds_ratio *= sample_value
                else:
                    odds_ratio /= sample_value
            sample_odds_ratios.append(log10(odds_ratio))
        sample_odds_ratios.sort()
        
        
    
        
        sample_odds_ratios_index = 0
        while sample_odds_ratios_index < number_simulations and sample_odds_ratios[sample_odds_ratios_index] <= 0:
            sample_odds_ratios_index += 1
        if sample_odds_ratios_index > number_simulations/2.0:
            sample_odds_ratios_index = number_simulations-sample_odds_ratios_index
        sample_odds_ratios_index *=2.0
        p_value = sample_odds_ratios_index / number_simulations
    
        
        edge_significant = False
        if p_value < alpha:
            edge_significant = True
            
        graph.add_edge(antibiotic_pair[0], antibiotic_pair[1], {'significant':edge_significant, 'odds_ratio':raw_odds_ratio})
    graphs_dict[source] = graph





for graph_name_pair in product(source_list, repeat=2):
    print graph_name_pair[0] + ':' + graph_name_pair[1]
    graph_2 = graphs_dict[graph_name_pair[1]]
    for edge in graphs_dict[graph_name_pair[0]].edges_iter(data=True):
        if edge[2]['significant']:
            graph_2_edge_data = graph_2.get_edge_data(edge[0], edge[1])
            if not graph_2_edge_data['significant']:
                print str(edge[0]) + ' ' + str(edge[1]) + ' ' + str(edge[2]['odds_ratio'])  + ' ' + str(graph_2_edge_data['odds_ratio']) 
        
    