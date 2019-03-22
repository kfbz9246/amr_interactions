'''
Created on Jul 5, 2016

@author: kfbz9246
'''

from numpy import array
from scipy.misc import comb


def zellens_test(table):
    for dimension in table.shape:
        if dimension !=2:
            print 'uh oh dimension not 2 test invalid'
    comb_1 = comb(sum(table[0,:,0]), table[0, 0, 0])
    comb_2 = comb(sum(table[0,:,1]), sum(table[0,0,:])-table[0,0,0])
    comb_3 = comb(sum(table[1,:,0]), sum(table[:,0, 0])-table[0,0,0])
    comb_4 = comb(sum(table[1,:,1]), sum(table[1,0,:])-sum(table[:,0,0])+table[0,0,0])
    
    print comb_1
    print comb_2
    print comb_3
    print comb_4
    
    numerator = comb_1*comb_2*comb_3*comb_4
    #print numerator

    
    
    z_min = max([0,sum(table[0,:,0])+sum(table[0,0,:])-sum(table[0].flatten()), sum(table[0,:,0])+sum(table[:,0,0])-sum(table[:,:,0].flatten())])
    z_max = min(sum(table[0,0]),sum(table[0,:,0]),sum(table[:,0,0]))+1
    print 'a'
    print z_min
    print z_max
    print 'b'
    
    denominator = 0
    for z in range(z_min, z_max):
        comb_1 = comb(sum(table[0,:,0]), z)
        comb_2 = comb(sum(table[0,:,1]), sum(table[0,0,:])-z)
        comb_3 = comb(sum(table[1,:,0]), sum(table[:,0, 0])-z)
        comb_4 = comb(sum(table[1,:,1]), sum(table[1,0,:])-sum(table[:,0,0])+z)
        denominator += comb_1*comb_2*comb_3*comb_4

        
    print denominator    
    return numerator/denominator

#def r_zellens_test(table):

table_1 = array([[100,200],[300,400]])
talbe_2 = array([[500,600],[700,800]])

table = array([table_1, talbe_2])

print table
print zellens_test(table)


