'''
Created on Dec 9, 2016

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
from networkx import Graph, draw_networkx_nodes, draw_networkx_labels, draw_networkx_edges, spring_layout, draw_networkx_edge_labels, draw_networkx, write_graphml
from itertools import combinations, izip, product
from numpy.random import dirichlet
import statsmodels.formula.api as smf
from statsmodels.api import families

antibiotic_set = ['cot', 'gen', 'str', 'fis', 'tet', 'kan', 'chl', 'amp'] #['cip', 'gen', 'fis']#['chl', 'tet', 'cip']#
dimention = 3
microbe_genus = 'EC'
host = 'Chickens'
start_year = 1995
end_year = 2013

print_pvalue = False
print_raw_data = False
print_contingency_table = False # table, row, column
print_raw_odds_ratio = False
print_loglinear_model = False
plot_posteriors = False
draw_graph = False
output_pvalue_spreadsheet = True

number_simulations = 1000000
alpha = .00001






engine = create_engine('sqlite:////home/code/amr_db/narms_database.sqlite')
session = sessionmaker(bind=engine)()

resistant_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('1').label('Resistance'), Samples.year.label('Year')).\
                         join(Isolates, Microbes, Samples).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic>=Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
suscptable_query = session.query(Tests.isolate_id.label('IsolateID'), Tests.drug_id.label('DrugID'), literal_column('0').label('Resistance'), Samples.year.label('Year')).\
                         join(Isolates, Microbes, Samples).\
                         join(Breakpoints,and_(Breakpoints.drug_id==Tests.drug_id, Breakpoints.genus==Microbes.genus)).\
                         filter(Tests.corrected_mic<Breakpoints.mic_r, Samples.year>=start_year, Samples.year<=end_year, Microbes.genus==microbe_genus, Samples.host_species==host)
    
if len(antibiotic_set) > 0:
    resistant_query = resistant_query.filter(Tests.drug_id.in_(antibiotic_set))
    suscptable_query = suscptable_query.filter(Tests.drug_id.in_(antibiotic_set))


data_query = resistant_query.union(suscptable_query)
data_all_years = read_sql(data_query.statement, engine)


all_pvalues_dict = {}
for year, data in data_all_years.groupby('Year'):
    pvalue_dict = {}

    data =  data.pivot('IsolateID', 'DrugID', 'Resistance')
    
    data.dropna(1,thresh=int(len(data)*.2), inplace=True)
    data.dropna(how='any',inplace=True)
    
    if len(antibiotic_set) >= 0:
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
    
    if draw_graph:
        graph = Graph()
        for antibiotic in antibiotic_set:
            graph.add_node(antibiotic)
    
    
    for antibiotic_pair in combinations(antibiotic_set, dimention):
        if print_pvalue or print_raw_data or print_contingency_table or print_raw_odds_ratio:
            print antibiotic_pair
                
        test_data = data[list(antibiotic_pair)+['Count']]
        test_data = test_data.groupby(antibiotic_pair).sum().reset_index()
        
        if print_raw_data:
            print test_data
    
        if print_contingency_table:
            data_as_array = zeros([2]*dimention,dtype=int)
            for index, row in test_data.iterrows():
                data_as_array[tuple(row.iloc[range(dimention)])] = int(row['Count'])
            print data_as_array
            
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
        
        if print_raw_odds_ratio:
            odds_ratio = 1
            for is_even, sample_value in zip(cell_counts_even, cell_counts_list):
                if is_even:
                    odds_ratio *= sample_value-.5
                else:
                    odds_ratio /= sample_value-.5
            print 'raw odds ratio:' + str(log10(odds_ratio))
        
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
        if print_pvalue:
            print 'P-value: '+ str(p_value)
        
        if output_pvalue_spreadsheet:
            pvalue_dict[antibiotic_pair] = p_value
        #print str(sample_odds_ratios[conficence_lower_bound]) + ' - ' + str(sample_odds_ratios[confidence_upper_bound])
        #print '--------------------------------'
        if draw_graph:
            if p_value < alpha:
                graph.add_edge(antibiotic_pair[0], antibiotic_pair[1])
        
        if plot_posteriors:
            temp_fig = figure(figsize=(9.25, 7), dpi=96, frameon=True)
            fig_ax = temp_fig.add_subplot(111)
            fig_ax.hist(sample_odds_ratios, bins=100)
            fig_ax.set_title('Posterior Distribution')
            fig_ax.set_xlabel('Difference in Log Odds Ratios')
            fig_ax.set_ylabel('Count')
            out_file_name = ''
            for antibiotic in antibiotic_pair:
                out_file_name += str(antibiotic) + '-'
            out_file_name += '.png'
            savefig(out_file_name, bbox_inches='tight')
          
        if print_loglinear_model:
            model_string = ''
            for column_name in antibiotic_pair:
                model_string += column_name +'+'
            model_string = model_string[0:-1]
            model = smf.glm(formula="Count ~ ("+ model_string +")**3", data=test_data, family=families.Poisson())
            model_fit = model.fit()
            print model_fit.summary()
    if output_pvalue_spreadsheet:
        all_pvalues_dict[year] = pvalue_dict
            
    if draw_graph:
        temp_fig = figure(figsize=(9.25, 7), dpi=96, frameon=True)
        fig_ax = temp_fig.add_subplot(111)
        draw_networkx(graph, ax=fig_ax)
        savefig('figure_1.png', bbox_inches='tight')
        write_graphml(graph,'graph.graphml')
if output_pvalue_spreadsheet:
    DataFrame(all_pvalues_dict).to_csv('coefs.csv')
    from pickle import dump
    dump(all_pvalues_dict, open('all_pvalues_dict.pickle', 'w'))