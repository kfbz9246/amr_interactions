'''
Created on Oct 7, 2016

@author: kfbz9246
'''
from sqlalchemy import create_engine, func, or_, and_
from sqlalchemy.sql.expression import literal_column
from database_interface.database_interface import Samples, Microbes, Isolates, Drugs, Breakpoints, Tests
from sqlalchemy.orm import sessionmaker, aliased
from pandas import read_sql, Series, read_csv
from numpy import mean, std, array, mean, var
from matplotlib.pyplot import figure, show, savefig, get_cmap
from math import sqrt, log, isnan
from numpy import sum
#import statsmodels.formula.api as smf
#from statsmodels.api import families

from itertools import combinations
from cPickle import dump, load
from itertools import izip

data = read_csv('AMR-long.csv')
data = data.set_index(['Year', 'ID', 'Drug',]).unstack().reset_index()
data.columns = [str(col[0]) if str(col[1]) == '' else str(col[1]) for col in data.columns.values]

data.insert(0,'Count', 1)
data = data[data['Year']==2011]
data = data.groupby(['GEN', 'KAN']).sum().reset_index()
print data
data.to_csv('will.csv')