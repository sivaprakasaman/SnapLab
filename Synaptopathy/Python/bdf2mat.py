# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:54:58 2019

@author: AndrewSivaprakasam
"""

import mne
import numpy as np # support for large, multi-dimensional arrays and metrices
from scipy.io import savemat

#Concatenation code from Rav
from EEGpp import EEGconcatenateFolder

######################################################################################################################################################
# Preprocessing
######################################################################################################################################################

#Consider making this a function at some point

########################################################################################################################################################


froot = "C:\\Users\\racqu\\Documents\\Research\\Purdue\\HumanData\\AS\\SQ25\\"
fs = 16384
    
topchans = [31] #CHANGE AS NEEDED

full_raw, full_eves = EEGconcatenateFolder(froot,nchans = 34, refchans = ['EXG1','EXG2'], exclude=[], fs_new = 4e3)


event_id = {'Positive':1}

epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = 0.0, tmax = 1.5, reject_tmax = 1.3, picks = topchans, reject = dict(eeg = 100e-6))
epochs.load_data()
epochs_filtered = epochs.filter(70,None)

pos_data = epochs.get_data()

event_id = {'Negative':2}

epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = 0.0, tmax = 1.5, reject_tmax = 1.3, picks = topchans, reject = dict(eeg = 100e-6))
epochs.load_data()
epochs_filtered = epochs.filter(70,None)

neg_data = epochs.get_data()

neg_l = len(neg_data)
pos_l = len(pos_data)

length = (neg_l>=pos_l)*pos_l + (neg_l<pos_l)*neg_l 

tot_arr = np.zeros(2* length, dtype = np.object)

ind = 0

for i in range(0,length):
    
    tot_arr[ind] = pos_data[i]
    tot_arr[ind+1] = neg_data[i]        
    ind = ind + 2

savemat('SQ25_Data_full.mat',{'SQ25_tot_arr':tot_arr})





