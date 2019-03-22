'''
Created on Aug 17, 2017

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

query_antibiotic_set = ['cot', 'gen', 'str', 'fis', 'tet', 'kan', 'chl', 'azi', 'cip', 'nal']#['amp', 'tio', 'amc', 'fox', 'axo']#['amp', 'tio', 'amc', 'fox', 'axo', 'chl']
dimention = 3
microbe_genus = 'EC'

source = None



start_year = 2011
end_year = 2013


number_simulations = 1000000
alpha = .01#.00001#






engine = create_engine('sqlite:////home/workspace/amr_db/narms_database.sqlite')


interaction_data = {}
for host in ['Chickens', 'Turkeys']:
    session = sessionmaker(bind=engine)()
    
    resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Drugs.name.label('DrugID'), literal_column('1').label('Resistance')).\
                         join(Isolates, Microbes, Samples, Drugs).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Drugs.id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
    suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Drugs.name.label('DrugID'), literal_column('0').label('Resistance')).\
                         join(Isolates, Microbes, Samples, Drugs).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Drugs.id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
    
    if len(query_antibiotic_set) > 0:
        resistant_query = resistant_query.filter(Drugs.id.in_(query_antibiotic_set))
        suscptable_query = suscptable_query.filter(Drugs.id.in_(query_antibiotic_set))
    
    if source != None:
        resistant_query = resistant_query.filter(Samples.source == source)
        suscptable_query = suscptable_query.filter(Samples.source == source)
        
        
    
    data_query = resistant_query.union(suscptable_query)
    data = read_sql(data_query.statement, engine)
    
     
    
    data =  data.pivot('IsolateID', 'DrugID', 'Resistance')
    data.dropna(1,thresh=int(len(data)*.2), inplace=True)
    data.dropna(how='any',inplace=True)
    

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
    
    
    
    odds_ratio_dict = {}
    pvalue_dict = {}
    for antibiotic_pair in combinations(antibiotic_set, dimention):
        interaction_str = ''
        for drug in antibiotic_pair:
            interaction_str += '-'.join((' '.join((y.capitalize() for y in x.split())) for x in drug.split('-'))) + ':'
        interaction_str = interaction_str[0:-1]
            
        print interaction_str
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
        
        odds_ratio = 1
        for is_even, sample_value in zip(cell_counts_even, cell_counts_list):
            if is_even:
                odds_ratio *= sample_value-.5
            else:
                odds_ratio /= sample_value-.5
        odds_ratio_dict[interaction_str]=log10(odds_ratio)
        
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
    
        
    
        pvalue_dict[interaction_str] = p_value
        
     
    interaction_data['odds_ratio_'+str(host)] = odds_ratio_dict
    interaction_data['p_value_'+str(host)] = pvalue_dict

interaction_data = DataFrame(interaction_data)
interaction_data =  interaction_data[(interaction_data['p_value_Chickens']<alpha) | (interaction_data['p_value_Turkeys']<alpha)]
interaction_data.to_csv(str(dimention)+'.csv')

    
    