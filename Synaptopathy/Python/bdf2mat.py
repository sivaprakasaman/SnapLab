# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:54:58 2019

@author: AndrewSivaprakasam
"""

import mne
import numpy as np # support for large, multi-dimensional arrays and metrices
from scipy.io import savemat
from anlffr.spectral import mtplv
from anlffr.spectral import mtspec
import matplotlib.pyplot as plt
import os

#Concatenation code from Rav
from EEGpp import EEGconcatenateFolder

######################################################################################################################################################
# Preprocessing
######################################################################################################################################################

#####################################################################################################################################################

def bdf2mat(froot,fs,Fs_new,hpf,t_stim,topchans,trial_name):
    full_raw, full_eves = EEGconcatenateFolder(froot,nchans = 34, refchans = ['EXG1','EXG2'], exclude=[], fs_new = Fs_new)
    
    event_id = {'Positive':1}
    
    epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = t_stim[0], tmax = t_stim[1], reject_tmax = 1.3, picks = topchans, reject = dict(eeg = 100e-6))
    epochs.load_data()
    epochs_filtered = epochs.filter(hpf,None)
    
    pos_data = epochs_filtered.get_data()
    
    event_id = {'Negative':2}
    
    epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = t_stim[0], tmax = t_stim[1], reject_tmax = 1.3, picks = topchans, reject = dict(eeg = 100e-6))
    epochs.load_data()
    epochs_filtered = epochs.filter(hpf,None)
    
    neg_data = epochs_filtered.get_data()
    
    neg_l = len(neg_data)
    pos_l = len(pos_data)
    
    length = (neg_l>=pos_l)*pos_l + (neg_l<pos_l)*neg_l 
    
    tot_arr = np.zeros(2* length, dtype = np.object)
    
    ind = 0
    
    for i in range(0,length):
        
        tot_arr[ind] = pos_data[i]
        tot_arr[ind+1] = neg_data[i]        
        ind = ind + 2
    
    os.chdir(froot)
    savemat(trial_name+'_Data_full.mat',{(trial_name+'_tot_arr'):tot_arr})


###########################PLOT SPECTROGRAM/PLV####################################

    x = np.add(neg_data[:length,:,:],pos_data[:length,:,:])
    
    params = dict(Fs = 4e3, tapers = [1,1],fpass = [0,1000], itc = 0)
    
    S_psd, N_psd, f_psd = mtspec(x,params)
    mt_plv, f_plv = mtplv(x,params)
    
    fig, ax = plt.subplots(num= 2, figsize=(12,8))
    ax.plot(f_psd, np.subtract(S_psd,N_psd))
    #ax.plot(f_psd, N_psd)
    ax.grid(color='k', linestyle='-', linewidth=1)
    ax.set_xlabel('Freq (Hz)', fontsize=18)
    ax.set_ylabel('PSD', fontsize=18)
    ax.set_xlim([70., 500.])
    ax.set_title('psd all')
    
    fig, ax = plt.subplots(num= 1, figsize=(12,8))
    ax.plot(f_plv, mt_plv)
    #ax.plot(f_psd, N_psd)
    ax.grid(color='k', linestyle='-', linewidth=1)
    ax.set_xlabel('Freq (Hz)', fontsize=18)
    ax.set_ylabel('PLV', fontsize=18)
    ax.set_xlim([70., 500.])
    ax.set_title('PLV all')
    
##############################################################################    
trial_name = 'SQ25'    
froot = "C:\\Users\\racqu\\Documents\\Research\\Purdue\\HumanData\\AS\\"+trial_name+'\\'
fs = 16384
Fs_new = 4e3
hpf = 70
t_stim = [0.0,1.5] 
topchans = [31] #CHANGE AS NEEDED

bdf2mat(froot,fs,Fs_new,hpf,t_stim,topchans,trial_name)




