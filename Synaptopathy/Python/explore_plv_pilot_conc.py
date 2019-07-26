#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:26:24 2019

@author: parida
"""

import mne
#from anlffr.helper import biosemi2mne as bs
from EEGpp import EEGconcatenateFolder
from anlffr.spectral import mtplv
from anlffr.spectral import mtspec
from glob import glob
import matplotlib.pyplot as plt
from os import mkdir
import os.path as op


dir_dict= {'dataIn':'\\Users\\racqu\\Documents\\Research\\Purdue\\HumanData\\AS\\SQ25\\', \
           'figOut':'\\Users\\racqu\\Documents\\Research\\Purdue\\HumanData\\Figs\\'} # define directories # define directories
read_data_params= {'nchans': 34, 'refchans':['EXG1','EXG2'], 'exclude': ['EXG3','EXG4','EXG5','EXG6','EXG7','EXG8']}

if not op.isdir(dir_dict['figOut']):
    mkdir(dir_dict['figOut'])
allFiles= glob(dir_dict['dataIn'] + '*.bdf') 
plotYes= True
saveYes= False


# Loop through directories
#    raw, eves= bs.importbdf(fileVar, refchans=read_data_params['refchans'], exclude=read_data_params['exclude'])
raw, eves= EEGconcatenateFolder(dir_dict['dataIn'], nchans=read_data_params['nchans'], fs_new=4e3, \
                                refchans=read_data_params['refchans'], exclude=read_data_params['exclude'])
#    raw= mne.io.read_raw_edf(fileVar)
#    raw.load_data()
#    eves= mne.find_events(raw, shortest_event=1, mask=255)
#    raw.filter(70, 1e3, phase='zero') # Not needed here as mtspec/mtplv have filter params
epochs= mne.Epochs(raw, eves, 1, tmin=-0.1, proj=False, tmax=1.2, baseline=(-0.1, 0.0), picks = 31, reject=dict(eeg=200e-6))
evoked= epochs.average()
params= dict(Fs=raw.info['sfreq'], tapers=[1, 1], fpass = [70, 1000], itc=0)
x = epochs.get_data()
x = x[:,:32, :]
x = x.transpose((1,0,2))
plv, f_plv = mtplv(x[30, :, :], params)
S_psd, N_psd, f_psd = mtspec(x[30, :, :], params)


if plotYes:
    fig,ax= plt.subplots(num=1, figsize=(12, 8))
    ax.plot(f_plv, plv)
    ax.grid(color='k', linestyle='-', linewidth=1)
    ax.set_xlabel('Freq (Hz)', fontsize=18)
    ax.set_ylabel('PLV', fontsize=18)
    ax.set_xlim([70., 500.])
    ax.set_title('plv all')
    if saveYes:
        plt.savefig(dir_dict['figOut']+'plv_all.png')
    
if plotYes:
    fig, ax = plt.subplots(num= 2, figsize=(12,8))
    ax.plot(f_psd, S_psd)
    ax.plot(f_psd, N_psd)
    ax.grid(color='k', linestyle='-', linewidth=1)
    ax.set_xlabel('Freq (Hz)', fontsize=18)
    ax.set_ylabel('PSD', fontsize=18)
    ax.set_xlim([70., 500.])
    ax.set_title('psd all')
    if saveYes:
        plt.savefig(dir_dict['figOut']+'psd_all.png')
