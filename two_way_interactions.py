'''
Created on Aug 24, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, merge, DataFrame
from numpy import mean, std, array, mean, var, zeros
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log10, isnan
from numpy import sum
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx
from itertools import combinations, izip, product
from numpy.random import dirichlet


plot_posteriors = False
conficence_lower_bound = 25
confidence_upper_bound = 999075
dimention = 3

antibiotic_set = ['chl', 'tet', 'cip']#['axo', 'tio', 'amp']#['chl', 'cip', 'gen', 'kan', 'fis', 'str', 'cot', 'tet']#['amc', 'amp', 'axo', 'chl', 'cot', 'fis', 'fox', 'gen','kan', 'nal', 'str', 'tet', 'tio']
microbe_genus = 'EC'

generate_new_graph_positions = False
show_order = False

engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=2011, Samples.year<=2013, Microbes.genus==microbe_genus, Samples.host_species=='Chickens')
suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance')).\
                     join(Isolates, Microbes, Samples).\
                     join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                     filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=2011, Samples.year<=2013, Microbes.genus==microbe_genus, Samples.host_species=='Chickens')

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

'''
graph = Graph()

for antibiotic in antibiotic_set:
    graph.add_node(antibiotic)
'''
test_results = {}
for antibiotic_pair in combinations(antibiotic_set, dimention):
    test_data = data[list(antibiotic_pair)+['Count']]
    test_data = test_data.groupby(antibiotic_pair).sum().reset_index()
    
    print test_data
    
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

    odds_ratio = 1
    for is_even, sample_value in zip(cell_counts_even, cell_counts_list):
        if is_even:
            odds_ratio *= sample_value-.5
        else:
            odds_ratio /= sample_value-.5
    print log10(odds_ratio)
    
    sample_odds_ratios = []
    for i in range(1000000):
        sample = dirichlet(cell_counts_list)
        odds_ratio = 1
        for is_even, sample_value in zip(cell_counts_even, sample):
            if is_even:
                odds_ratio *= sample_value
            else:
                odds_ratio /= sample_value
        sample_odds_ratios.append(log10(odds_ratio))

    sample_odds_ratios.sort()
    
    test_results[antibiotic_pair] = {'2.5_quantile':sample_odds_ratios[conficence_lower_bound], '97.5_quantile':sample_odds_ratios[confidence_upper_bound], 'significance':(0 if sample_odds_ratios[conficence_lower_bound] < 0 and 0< sample_odds_ratios[confidence_upper_bound] else 1)}
    '''
    if not(sample_odds_ratios[conficence_lower_bound] < 0 and 0< sample_odds_ratios[confidence_upper_bound]):
        graph.add_edge(antibiotic_pair[0], antibiotic_pair[1])
    '''
    print antibiotic_pair
    sample_odds_ratios_index = 0
    while sample_odds_ratios[sample_odds_ratios_index] <= 0:
        sample_odds_ratios_index += 1
    if sample_odds_ratios_index > 500000:
        sample_odds_ratios_index = 1000000-sample_odds_ratios_index
    sample_odds_ratios_index *=2
    sample_odds_ratios_index /= 1000000.0
    print 'P-value: '+ str(sample_odds_ratios_index)
    print str(sample_odds_ratios[conficence_lower_bound]) + ' - ' + str(sample_odds_ratios[confidence_upper_bound])
    print '--------------------------------'
    if plot_posteriors:
        temp_fig = figure(figsize=(9.25, 7), dpi=96, frameon=True)
        fig_ax = temp_fig.add_subplot(111)
        fig_ax.hist(sample_odds_ratios, bins=100)
        fig_ax.set_title('Posterior Distribution Of Ratio Of Odds')
        out_file_name = ''
        for antibiotic in antibiotic_pair:
            out_file_name += str(antibiotic) + '-'
        out_file_name += '.png'
        savefig(out_file_name, bbox_inches='tight')
DataFrame(test_results).to_csv('test_results.csv')
'''
temp_fig = figure(figsize=(9.25, 7), dpi=96, frameon=True)
fig_ax = temp_fig.add_subplot(111)
draw_networkx(graph, ax=fig_ax)
savefig('figure_1.png', bbox_inches='tight')

from pickle import dump
dump(graph, open('graph.pickle', 'w'))
'''
