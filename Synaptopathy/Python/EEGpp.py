# -*- coding: utf-8 -*-
"""
Created on Sat Sep 08 21:30:13 2018

@author: StuffDeveloping
EEG pre-processing
"""
from os import listdir
from anlffr.helper import biosemi2mne as bs
from mne import concatenate_raws
from matlab2python import numel
import numpy as np


def EEGconcatenateFolder(folder,nchans, refchans, fs_new=[], exclude=[]):
    #concatenates all the EEG files in one folder
    #Assumes files are in format output from biosemi ... subj.bdf, subj+001.bdf, etc.
    #Also folder should end with a '/' 
    EEGfiles = listdir(folder)
    EEGfiles.sort()# This line and next to fix order of files 
    EEGfiles.insert(0,EEGfiles.pop(len(EEGfiles)-1))
    print(EEGfiles)
    raw = [] 
    events = []
    
    for eeg_f in EEGfiles:
        raw_temp, events_temp= bs.importbdf(folder+eeg_f,nchans,refchans,exclude =exclude) #uses EXG1 and EXG2 as reference usually
        
        if numel(fs_new):
            print('Resample raw data and update eves indices')
            events_temp[:,0]= np.round(events_temp[:,0]/raw_temp.info['sfreq']*fs_new).astype('int')
            raw_temp.resample(fs_new)
            
        raw.append(raw_temp)
        events.append(events_temp)
    EEG_full, events_full = concatenate_raws(raw,events_list=events)
    return EEG_full, events_full