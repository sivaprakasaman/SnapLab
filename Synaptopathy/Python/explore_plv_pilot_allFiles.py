#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:26:24 2019

@author: parida
"""

import mne
from anlffr.helper import biosemi2mne as bs
#from EEGpp import EEGconcatenateFolder
from anlffr.spectral import mtplv
from glob import glob
import matplotlib.pyplot as plt
from os import mkdir
import os.path as op


dir_dict= {'dataIn':'\Users\racqu\Documents\Research\Purdue\Human Data\\', \
           'figOut':'\Users\racqu\Documents\Research\Purdue\Human Data\Figs\\'} # define directories
read_data_params= {'refchans':['EXG1','EXG2'], 'exclude': ['EXG3','EXG4','EXG5','EXG6','EXG7','EXG8']}

if not op.isdir(dir_dict['figOut']):
    mkdir(dir_dict['figOut'])
allFiles= glob(dir_dict['dataIn'] + '*.bdf') 
plotYes= True


# Loop through directories
for fileVar in allFiles:
    print('Working on ->' + fileVar)
    raw, eves= bs.importbdf(fileVar, refchans=read_data_params['refchans'], exclude=read_data_params['exclude'])
#    raw= mne.io.read_raw_edf(fileVar)
#    raw.load_data()
#    eves= mne.find_events(raw, shortest_event=1, mask=255)
    raw.filter(70, 1e3, phase='zero')
    epochs= mne.Epochs(raw, eves, 1, tmin=-0.1, proj=False,tmax=1.2, baseline=(-0.1, 0.0),reject=dict(eeg=200e-6))
    evoked= epochs.average()
    params= dict(Fs=raw.info['sfreq'], tapers=[1, 1], fpass = [70, 1000], itc=0)
    x = epochs.get_data()
    x = x[:,:32, :]
    x = x.transpose((1,0,2))
    plv, f = mtplv(x[30, :, :], params)
    
    
    if plotYes:
        filename= op.basename(fileVar)
        filename, _= op.splitext(filename)
        fig,ax= plt.subplots(num=1, figsize=(12, 8))
        ax.plot(f, plv)
        ax.grid(color='k', linestyle='-', linewidth=1)
        ax.set_xlabel('Freq (Hz)', fontsize=18)
        ax.set_ylabel('PLV', fontsize=18)
        ax.set_xlim([70., 500.])
        ax.set_title(filename)
        plt.savefig(dir_dict['figOut']+filename+'.png')