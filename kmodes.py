import pandas as pd
import sys
import random
import math
import re
import matplotlib.pyplot as plt
import numpy as np
import os
import time


#initialize centers
def initialize(data, cols):
    print('The initilizing clusters are:')
    init = []
    i=0
    while(i < km):
        u = []
        for j in range(len(cols)):
            if cols[j] == 'CLNSIG':
                chc = data[cols[j]].tolist()
                chc = [x for x in chc if (x != 'Uncertain_significance' and x != 'not_provided' )]
                u.append(random.choice(chc))
            else:
                u.append(random.choice(data[cols[j]].tolist()))
        if u not in init:
            init.append(u)
            i+=1
            print(u)
    return init

#update clusters-updates location of the means to the centroid of each cluster
def upcl(modes, data):
    for md_idx in range(len(modes)): 
        dataT = data.loc[data['clstrs'] == md_idx]
        if len(dataT) != 0:
            new_mode = dataT.mode()
            mode = new_mode.iloc[0][:-1].values.tolist().copy()
            if mode not in modes:
                modes[md_idx] = mode
    return modes

#adds the total distance from all features of a point Uncertain_significance
def vectD(x,y,cols):
    count = 0
    for i in range(len(x)):
        if y[i] != 'Uncertain_significance':
            if(x[i]==y[i]):
                if cols[i] == 'CLNDN' or 'CLNVC':
                    count += 2
                elif cols[i] == 'CLNSIG':
                    count += 3
                else:
                    count += 1
    return count

#check distance function checks the closest cluster to each point and 
#   returns a closest-cluster array
def chk_D(modes, data):
    clstrs = [None]*len(data)
    for j in range(len(data)):
        vect = data.iloc[[j]].values.tolist()
        vect = vect[0][:]
        D= [None]*len(modes) 
        for k in range(len(modes)):
            D[k] = vectD(modes[k],vect,data.columns)
            if(k==0):
                cl = k
            elif(D[cl]>D[k]):
                cl = k
        clstrs[j]=cl
    return clstrs

#Apply k-means algorithm
def kmodes(data,cols,t,max_iter):
    #step 1 initialize
    N, step, conv = .000001,0,1
    init = initialize(data, cols)
    data['clstrs'] = ""
    modes = init.copy()
    #iterate step 2 and 3 until convergence or 100 iterations
    while((conv>N) and (step <max_iter)):
        #step 2
        data['clstrs'] = chk_D(modes, data.drop(columns=['clstrs']))
        #step 3
        if(step>1):
            old_modes2 = old_modes.copy()
        else:
            old_modes2 = modes.copy()
        old_modes = modes.copy()
        modes = upcl(modes, data)
        #convergence condition is if they are equal or if they are alternating between 2 mode sets
        if((sorted(modes)==sorted(old_modes)) or (sorted(modes)==sorted(old_modes2))):
            break
        step+=1
        elapsed = time.time() - t
        print('Iteration',step,'completed in',round(elapsed/60), 'minutes and', round(elapsed%60), 'seconds.')
            
    return modes, data, step, init

def printer(modes, data,elapsed, iterations):
    print('Number of iterations ran is', iterations)
    print('The centroids are:')
    for i in range(len(modes)):
        num = len(data.loc[data['clstrs'] == i])
        if num >0:
            print('Centroid', i+1, 'is:',modes[i],'with', num, 'samples assigned to it.')
    print('The elapsed time is:',round(elapsed/60), 'minutes and', round(elapsed%60), 'seconds.')
    print('done.')

########################
########  MAIN  ########
########################

#prevariables, with default on my computer, this took ~15 min and converged at 196 iterations
os.chdir("E:\\Users\\Nathanael\\Documents\\")
filenm = 'Clinvardata_formated.tsv'
km = 10
max_iter = 500
samp_size = 1000

#function
#read in the raw data from file
rawDat = pd.read_csv(filenm, sep='\t', header = 0, dtype={'REF':'category','ALT':'category','CLNDN':'category','CLNVC':'category','MC':'category','GENEINFO':'category'})
rawcols = ['CHROM', 'POS','ID','REF','ALT','CLNSIG','CLNDN','CLNVC','MC','GENEINFO']
rawDat = rawDat[rawcols]
#sample 500 random data
data = rawDat.sample(n=samp_size)
data = data.drop(columns = ['CHROM','ID','POS'])
cols = data.columns.tolist()
t = time.time()
#run kmodes
modes, data, iterations, init = kmodes(data,cols, t, max_iter)
elapsed = time.time() - t
printer(modes,data,elapsed,iterations)

