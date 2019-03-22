'''
Created on Nov 4, 2016

@author: kfbz9246
'''
from pandas import DataFrame
from math import log10
from numpy import zeros
from numpy.random import dirichlet
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, savefig, get_cmap
import statsmodels.formula.api as smf
from statsmodels.api import families

dimention = 3
print_pvalue = True
print_raw_data = True
print_contingency_table = True # table, row, column
print_raw_odds_ratio = True
plot_posteriors = True
print_loglinear_model = False

number_simulations = 1000000
alpha = .01

antibiotic_pair = ('Z', 'X', 'Y',)

'''
test_data = [{'X':0, 'Y':0, 'Count':190},
             {'X':0, 'Y':1, 'Count':100},
             {'X':1, 'Y':0, 'Count':10},
             {'X':1, 'Y':1, 'Count':100}]
'''

test_data = [{'Z':0,'X':0, 'Y':0, 'Count':1},
             {'Z':0,'X':0, 'Y':1, 'Count':9},
             {'Z':0,'X':1, 'Y':0, 'Count':1},
             {'Z':0,'X':1, 'Y':1, 'Count':1},
             {'Z':1,'X':0, 'Y':0, 'Count':1},
             {'Z':1,'X':0, 'Y':1, 'Count':1},
             {'Z':1,'X':1, 'Y':0, 'Count':9},
             {'Z':1,'X':1, 'Y':1, 'Count':1}]






test_data = DataFrame(test_data)
test_data = test_data[['Z','X', 'Y', 'Count']]

if print_raw_data:
    print test_data
    test_data.to_csv('test_data.csv')

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


if plot_posteriors:
    temp_fig = figure(figsize=(9.25, 4.75), dpi=96, frameon=True)#(9.25, 7)
    fig_ax = temp_fig.add_subplot(111)
    fig_ax.hist(sample_odds_ratios, bins=100)#, color='.75')
    fig_ax.set_title('Posterior Distribution', fontsize=18)
    fig_ax.set_xlabel('Difference in Log Odds Ratios', fontsize=14)
    fig_ax.set_ylabel('Count', fontsize=14)
    fig_ax.tick_params(axis='both', which='major', labelsize=12)
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