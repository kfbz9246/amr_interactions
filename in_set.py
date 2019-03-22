'''
Created on Oct 20, 2016

@author: kfbz9246
'''
set_items = ['gen', 'str', 'cot', 'kan', 'fis', 'tet', 'amp', 'chl']#['amp', 'chl', 'axo', 'fox', 'amc', 'tio']#
the_set = {}
for item in set_items:
    the_set[item] = None

for row in open('foo.txt'):
    row_in_set = True
    for item in row.split():
        if not the_set.has_key(item):
            row_in_set = False
    if row_in_set:
        print row
    