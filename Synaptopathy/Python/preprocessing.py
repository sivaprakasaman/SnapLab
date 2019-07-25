# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:54:58 2019

@author: AndrewSivaprakasam
"""

from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np # support for large, multi-dimensional arrays and metrices
import os # it assigns its path attribute to an os-specific path module
from scipy.io import savema
from EEGpp import EEGconcatenateFolder
 
#import xlwt
#from mne.time_frequency import tfr_multitaper
#from mne.connectivity import spectral_connectivity
#import fnmatch # unix filename pattern matching
#from scipy.signal import butter, filtfilt, hilbert  
#import pylab as pl
#from anlffr.spectral import mtplv


def offsetSub(evoked, fs):
    offset = evoked[int(0.525*fs):int(0.555*fs)+1] # add 0.025 s
    subtraction = np.concatenate((np.zeros(int(0.275*fs)), offset)) # add 0.025 s 
    subtraction =  np.concatenate((subtraction, np.zeros(len(evoked)-len(subtraction))))
    return subtraction 

def highpass_filter(data, cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    y = filtfilt(b, a, data)
    return y
######################################################################################################################################################
# Preprocessing
######################################################################################################################################################
def evoked(raw, eves, eid):
    eegReject = 50e-6
    epochs = mne.Epochs(raw, eves, eid, tmin = 0.0016, proj = False, tmax = 0.0016 + 0.351, 
                                baseline = (None, None), reject = dict(eeg=eegReject)) # baseline correction is not that necessary since the data was already high-passed
    epochs_prob = mne.Epochs(raw, eves, eid, tmin = 0.0016, proj = False, tmax = 0.1 + 0.0016, 
                                baseline = (None, None), reject = dict(eeg=eegReject))
    epochs_masked = mne.Epochs(raw, eves, eid, tmin = 0.251 + 0.0016, proj = False, tmax = 0.0016 + 0.35105, 
                                baseline = (None, None), reject = dict(eeg=eegReject)) # picked 0.35105 instead of 0.351 to equalized the number of points of epochs_prob and epochs_masked 
    epochs_masker = mne.Epochs(raw, eves, eid, tmin = 0.15 + 0.0016, proj = False, tmax = 0.0016 + 0.25, 
                                baseline = (None, None), reject = dict(eeg=eegReject))
    data_prob = epochs_prob.get_data()
    data_masked = epochs_masked.get_data()
    data_masker = epochs_masker.get_data()
    
    numTrial = len(epochs.events) # number of trials
    evoked = epochs.average()
#    evoked.plot(picks=[30, 31])
    #topchans = [3, 4, 7, 22, 26, 25, 30, 31]
    tiptrodes = [32, 33]
    chans = topchans
    #evoked_all = evoked.data[chans, :].mean(axis = 0) - evoked.data[tiptrodes, :].mean(axis=0)
    evoked_all = evoked.data[chans, :].mean(axis = 0)
    return evoked_all, numTrial, data_prob, data_masked, data_masker
 
def processing(fpath, bdf, eid, fs): # eid = [id_475_pos, id_475_neg, id_500_pos, id_500_neg, id_525_pos, id_525_neg,]
    rawtemp, evestemp = bs.importbdf(fpath + bdf)
    raw = rawtemp
#    raw.info['bads'] += ['A24', 'A16', 'A25']  
    evestemp[:, 1] = np.mod(evestemp[:, 1], 256) 
    evestemp[:, 2] = np.mod(evestemp[:, 2], 256)
    raw.filter(l_freq = 400, h_freq=1300) # adjust this range 
    evoked_pos = []
    evoked_neg = []
    
    numpos = []
    numneg = []
    
    flags = [False, False]
    if eid[0] in evestemp[:,2]:
        evoked_pos, numpos, epochs_prob_pos, epochs_masked_pos, epochs_masker_pos = evoked(raw, evestemp, eid[0])
#        subtraction = offsetSub(evoked_pos, fs)
#        evoked_pos = evoked_pos - subtraction
#        evoked_pos = evoked_pos[int(0.0266*fs):int((0.0266+0.395)*fs)] # starting from time 0 to the end of masked stimulus   
        flags[0] = True;
        length = len(evoked_pos)
    if eid[1] in evestemp[:,2]:
        evoked_neg, numneg, epochs_prob_neg, epochs_masked_neg, epochs_masker_neg = evoked(raw, evestemp, eid[1])
#        subtraction = offsetSub(evoked_neg, fs)
#        evoked_neg = evoked_neg - subtraction 
#        evoked_neg = evoked_neg[int(0.0266*fs):int((0.0266+0.395)*fs)]   
        flags[1] = True;
        length = len(evoked_neg)
    
    return evoked_pos, evoked_neg, numpos, numneg, length, flags, epochs_prob_pos, epochs_masked_pos, epochs_masker_pos, epochs_prob_neg, epochs_masked_neg, epochs_masker_neg    

def weightedAvg(evoked, numTrial):
    numTotal = np.sum(numTrial)
    evokedAvg = np.zeros(len(evoked[1]))
    for k in range(len(numTrial)):
        evokedAvg = evokedAvg + evoked[k]*numTrial[k]/numTotal
    return evokedAvg

def evokedExtraction(froot, subject, trigNum, fs):
    fpath = froot + '/' + subject + '/'
    print('Running subject', subject)
    
    # extracting bdf filenames into bdfs
    bdfs = fnmatch.filter(os.listdir(fpath), subject + '_FWMK*.bdf') 

    numTrialpos = np.zeros(len(bdfs))
    numTrialneg = np.zeros(len(bdfs))
    
    if len(bdfs) >= 1:
        for k, bdf in enumerate(bdfs):                
            pos, neg, numpos, numneg, length, flags, epochs_prob_pos, epochs_masked_pos, epochs_masker_pos, epochs_prob_neg, epochs_masked_neg, epochs_masker_neg = processing(fpath, bdf, trigNum, fs) #trigNum: [1, 2], for example
            if k == 0:
                evoked_poss = np.zeros((len(bdfs), length))
                evoked_negs = np.zeros((len(bdfs), length))
                data_prob_pos = epochs_prob_pos
                data_prob_neg = epochs_prob_neg
                data_masked_pos = epochs_masked_pos
                data_masked_neg = epochs_masked_neg
                data_masker_pos = epochs_masker_pos
                data_masker_neg = epochs_masker_neg
            if k > 0:
                data_prob_pos = np.concatenate((epochs_prob_pos, data_prob_pos), axis = 0)
                data_prob_neg = np.concatenate((epochs_prob_neg, data_prob_neg), axis = 0)
                data_masked_pos = np.concatenate((epochs_masked_pos, data_masked_pos), axis = 0)
                data_masked_neg = np.concatenate((epochs_masked_neg, data_masked_neg), axis = 0)
                data_masker_pos = np.concatenate((epochs_masker_pos, data_masker_pos), axis = 0)
                data_masker_neg = np.concatenate((epochs_masker_neg, data_masker_neg), axis = 0)
            if flags[0]:    
                evoked_poss[k] = pos
                numTrialpos[k] = numpos
            if flags[1]:
                evoked_negs[k] = neg
                numTrialneg[k] = numneg
        
        evoked_pos = weightedAvg(evoked_poss, numTrialpos) 
        evoked_neg = weightedAvg(evoked_negs, numTrialneg) 
    else:
        RuntimeError("No bdf files found!")
    return evoked_pos, evoked_neg, data_prob_pos, data_prob_neg, data_masked_pos, data_masked_neg, data_masker_pos, data_masker_neg

def diff(epochs_pos, epochs_neg):
    numEvntsPos = epochs_pos.shape
    numEvntsNeg = epochs_neg.shape
    if numEvntsPos[0] > numEvntsNeg[0]:
        subt = epochs_pos[0:numEvntsNeg[0], :, :] - epochs_neg
    else:
        subt = epochs_pos - epochs_neg[0:numEvntsPos[0], :, :]
    return subt

def summ(epochs_pos, epochs_neg):
    numEvntsPos = epochs_pos.shape
    numEvntsNeg = epochs_neg.shape
    if numEvntsPos[0] > numEvntsNeg[0]:
        summation = epochs_pos[0:numEvntsNeg[0], :, :] + epochs_neg
    else:
        summation = epochs_pos + epochs_neg[0:numEvntsPos[0], :, :]
    return summation, min(numEvntsPos[0], numEvntsNeg[0])

def diffAdption(data_prob_pos, data_prob_neg, data_masked_pos, data_masked_neg):
    data_prob_diff = diff(data_prob_pos, data_prob_neg)
    data_masked_diff = diff(data_masked_pos, data_masked_neg)
    data_adpt_diff = diff(data_prob_diff, data_masked_diff)
    numtrials_adpt = data_adpt_diff.shape
    return data_adpt_diff, numtrials_adpt[0]

def subjectProcessing(froot, subj, eid, fs):
    evoked_pos, evoked_neg, data_prob_pos, data_prob_neg, data_masked_pos, data_masked_neg, data_masker_pos, data_masker_neg = evokedExtraction(froot, subj, eid, fs)
    data_adpt_diff, numtrials_adpt = diffAdption(data_prob_pos, data_prob_neg, data_masked_pos, data_masked_neg)
    # saving epochs and evoked
    np.savez(froot+'/'+'epochs_evoked_FCz' + '/'+subj, evoked_pos = evoked_pos, evoked_neg = evoked_neg, data_prob_pos = data_prob_pos, data_prob_neg = data_prob_neg, data_masked_pos = data_masked_pos, data_masked_neg = data_masked_neg, data_masker_pos = data_masker_pos, data_masker_neg = data_masker_neg, data_adpt_diff = data_adpt_diff) 
    data_prob_diff = diff(data_prob_pos, data_prob_neg)
    data_masker_diff = diff(data_masker_pos, data_masker_neg)
    data_sumMasker, _ = summ(data_masker_pos, data_masker_neg)
    return evoked_pos, evoked_neg, data_sumMasker, data_masker_diff, data_prob_diff, data_adpt_diff, numtrials_adpt
########################################################################################################################################################


froot = "C:\\Users\\racqu\\Documents\\Research\\Purdue\\HumanData\\AS\\SQ50\\"
fs = 16384
    
topchans = [31] #CHANGE AS NEEDED

bdfname = froot+bdf

full_raw, full_eves = EEGconcatenateFolder(froot,nchans = 34, refchans = ['EXG1','EXG2'], exclude=[], fs_new = 4e3)


event_id = {'Positive':1}

epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = 0.0, tmax = 1.5, reject_tmax = 1.3, picks = 31, reject = dict(eeg = 100e-6))
epochs.load_data()
epochs_filtered = epochs.filter(70,None)

pos_data = epochs.get_data()

event_id = {'Negative':2}

epochs = mne.Epochs(full_raw, full_eves, event_id,tmin = 0.0, tmax = 1.5, reject_tmax = 1.3, picks = 31, reject = dict(eeg = 100e-6))
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

savemat('SQ50_Data_full.mat',{'SQ50_tot_arr':tot_arr})





