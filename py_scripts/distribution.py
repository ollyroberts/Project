#!/usr/bin/python
import os
import sys
import re
import numpy as np
import scipy.stats as stats 
from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import mannwhitneyu

import matplotlib.pyplot as plt

import argparse



# generate gaussian data
from numpy.random import seed
from numpy.random import randn
from numpy import mean
from numpy import std
from numpy import array
from numpy import median


"""
The purpose of this program is to extract the xyz coordinates of the first (a)
 middle proline(b) and last (c) ca atom of protein

 if there is more than one proline at same dist from mid point select the lower number
"""

def main ():
    win_or_linux()

def win_or_linux():
    """This detects whether the system is windows or linux and then calls
    either windows_arguments() or linux_arguments() which controll how
    a file is opened """

    if sys.platform =='win32':
        file = windows_arguments()

    if sys.platform =='linux2':
        file = linux_arguments()

def windows_arguments():
    master('one_helix_angle.txt','two_helix_angle.txt','distribution_output.txt')

		

def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('one_helix', help = 'ca atoms of proteins seperated by chains')
    parser.add_argument('two_helix', help = 'output file which will have the pdb name and 1st residue with ABC angle')

    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    one_helix_file = vars(args)['one_helix']
    two_helix_file = vars(args)['two_helix']

    master(one_helix_file,two_helix_file)


def fileread(filename):
    """ opens a file object, removes the contents and strips lhwhitespace
    """

    file = open(filename, 'r')
    output = file.read()
    output = output.lstrip()
    file.close()

    return(output)

def angle_ext(input_file):

    angle = []

    helix_file = fileread(input_file)
    helix_file = helix_file.splitlines()

    counter = 0 
    for line in helix_file:
        length =(len(line.split()))
        split_line = line.split()
        angle.append((split_line[length - 1]))

    return angle

    #print(one_helix_angle)

def is_normal(input_list): 
    """
    input : list of strings e.g ['161.30411111', '148.4969377']
    """


    # converts from a list of str to a list of float 
    input_float = [float(i) for i in input_list]
    for x in input_float:
        print (x)

    # converts from a list to a numpy.ndarray
    data = array(input_float)


    mean_data       = mean(data)
    stnd_data       = std(data)
    skew_data       = skew(data)
    kurtosis_data   = kurtosis(data)

    # skewness is if the data is more on one side of the mean than the other
    # kurtosis score of 3 means the tails are normally waiting for a mean distribution.
    # higher indicates more weight in the tails.
    print('mean=' + str(mean_data)+' stdv=' +str(stnd_data) + ' skew=' +str(skew_data)+'kurtosis=' +str(kurtosis_data))
    print('(skewteset) s*2 + k*2 (kurtosistest), a2 sided chi squard probability for hypthesis :' + str(stats.normaltest(data)))
    #plt.style.use('ggplot')
    #plt.hist(data, bins=60)
    plt.title('Two helix bend angles with probability function')
    plt.xlabel('Angle of 6H(1-3Gap)6H')
    plt.ylabel('Frequency')

    # draw the probability function over it
    s = data 
    mu = mean_data
    sigma = stnd_data

    # generate univariate observations
    # summarize
    count, bins, ignored = plt.hist(s, 72,  density=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')

    #prints graph
    plt.show()







def manwhitneyu(list1, list2):
    """
    input: two lists 
    This is functing for the mann whitney u test which checks to see if tho
     sets have the same distribution
    """

    # converts from a list of str to a list of float

    list1_float = [float(i) for i in list1]
    list2_float = [float(i) for i in list2]

    # converts from a list to a numpy.ndarray
    data1 = array(list1_float)
    print(median(list1_float))
    data2 = array(list2_float)
    print(median(list2_float))

    print(mannwhitneyu(data1,data2))

#def skewness_and_kurtosis (data):

def qq_plot(list1):
    list1_float = [float(i) for i in list1]
    data1 = array(list1_float)

    N = 1200
    p = 0.53
    q = 1000
    obs = data1



    z = (obs-np.mean(obs))/np.std(obs)

    stats.probplot(z, dist="norm", plot=plt)
    plt.title("Normal Q-Q plot for two helix (6H1-3G6H) data")
    plt.show()

def t_test(list1, list2):
    """
    input: two lists 
    This is functing for the t-test which check if there is a correlationbetween
    two groups 
    """
    list1_float = [float(i) for i in list1]
    list2_float = [float(i) for i in list2]

    # converts from a list to a numpy.ndarray
    data1 = array(list1_float)
    data2 = array(list2_float)

    print(stats.ttest_ind(data1,data2))




def master(one_helix_input,two_helix_input,output_file):



    # unpacks and extracts the bond angles
    one_helix_list = angle_ext(one_helix_input)
    two_helix_list = angle_ext(two_helix_input)


    #checks to see if data is normally distributed 
    #is_normal(one_helix_list)
    #is_normal(two_helix_list)
    
    #manwhitneyu(one_helix_list,two_helix_list)

    #qq_plot(two_helix_list)
    t_test(one_helix_list,two_helix_list)




if __name__== "__main__":
    main()


