#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 16:25:08 2019

@author: agudemu
"""

from anlffr.helper import biosemi2mne as bs
from anlffr.spectral import mtplv
import mne
import numpy as np # support for large, multi-dimensional arrays and metrices
import pylab as pl
import os # it assigns its path attribute to an os-specific path module
import fnmatch # unix filename pattern matching
from scipy.signal import butter, filtfilt, hilbert  
from scipy.io import savemat  
from mne.time_frequency import tfr_multitaper
from mne.connectivity import spectral_connectivity 
import xlwt

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
    print 'Running subject', subject
    
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
# Data analysis
########################################################################################################################################################
def diffSumEvd(evoked_pos, evoked_neg):
    # evoked response analysis
    diffWhole = evoked_pos - evoked_neg
    PosProb = evoked_pos[0:int(0.1*fs)]
    NegProb = evoked_neg[0:int(0.1*fs)]
    diffProb = PosProb - NegProb
    sumProb = PosProb + NegProb
    PosMasked = evoked_pos[int(0.251*fs):int(0.351*fs)]
    NegMasked = evoked_neg[int(0.251*fs):int(0.351*fs)]
    diffMasked = PosMasked - NegMasked
    sumMasked = PosMasked + NegMasked
    PosMasker = evoked_pos[int(0.15*fs):int(0.25*fs)]
    NegMasker = evoked_neg[int(0.15*fs):int(0.25*fs)]
    diffMasker = PosMasker - NegMasker
    sumMasker = PosMasker + NegMasker
    adptPos = PosProb - PosMasked
    adptNeg = NegProb - NegMasked
    diffAdpt = adptPos - adptNeg
    sumAdpt = adptPos + adptNeg
    return diffWhole, diffProb, sumProb, diffMasked, sumMasked, diffMasker, sumMasker, diffAdpt, sumAdpt

def freqAnsis(sig):
    sig_fft = np.fft.fft(sig)
    magSig = np.abs(sig_fft)
    phase = np.angle(sig_fft)
    return magSig, phase

def fftPlot(evoked, fs, chop, chop2, stimType):
    # frequency analysis
    t = np.arange(0, len(evoked)/float(fs), 1/float(fs))
    mag, phase500 = freqAnsis(evoked[np.logical_and(t > chop, t < chop2)])
    y = evoked[np.logical_and(t > chop, t < chop2)]
    x = t[np.logical_and(t > chop, t < chop2)]
#    y = evoked
    
    freq = np.linspace(0, fs/4, len(mag)/4)
    fig = pl.figure(figsize = (20, 5))
    ax = fig.add_subplot(211)
#    pl.plot(t, diffWhole, label = 'prob')
    titleStr = subj + '_' + stimType
    pl.title(titleStr, fontsize=14)
    pl.xlabel('Time (s)', fontsize=14)
    pl.ylabel('Evoked response', fontsize=14)
    ax.plot(x, y)
    ax = fig.add_subplot(212)
    peak = max(mag[0:len(mag)/4])
    index = np.where(mag[0:len(mag)/4] == peak)
    ax.plot(freq, mag[0:len(mag)/4])
    ax.plot(freq[index], peak, 'r+', markersize = 12, linewidth = 8)
    ax.annotate(repr(peak), xy = (freq[index], peak), xytext = (freq[index], peak))
    pl.xlabel('Freq (Hz)', fontsize=14)
    pl.ylabel('Amplitude', fontsize=14)
    pl.savefig(froot + '/FWMKfigures_FCz/' + titleStr + '.png')
#    # saving mags
#    dictMat = {"mag": peak}
#    savemat(froot+'/mags/'+subj+stimType, dictMat)
    return peak, index

def plvAnalysis(epochs, Type):
    # phase locking analysis
    params = dict(Fs = 16384, tapers = [1, 1], fpass = [400, 600], itc = 0)
    # plv to TFS
    epochs = np.transpose(epochs, (1, 0, 2)) # switching the first and second columns
    plv, f = mtplv(epochs, params)
    index = np.where(plv[topchans[0],] == max(plv[topchans[0],])) 
    plv_max = max(plv[topchans[0],])
    f_max = f[index]
    dimen = epochs.shape
    numtrials = dimen[1]
    # saving plvs
    np.savez(froot+'/plvs/'+subj+Type, plv = plv, plv_max = plv_max, f = f_max, numtrials = numtrials);
    # collecting individual plvs into array for all subjects
    return plv_max, f_max, numtrials
    
########################################################################################################################################################
OS = 'Ubuntu'

if OS == 'Ubuntu':
    froot = '/media/agudemu/Storage/Data/EEG/FFR'
else:
    froot = '/Users/baoagudemu1/Desktop/Lab/EEG-Python/FFR'
#subjectList = ['S025', 'S031', 'S043', 'S051', 'S072', 'S075', 'S084', 'S117', 'S127', 'S128', 'S132', 'S133', 'S149', 'S183', 'S185', 
#               'S187', 'S191', 'S194', 'S195', 'S196', 'S197', 'S216'] # exlcuded S078, S123 (empty arrays were returned), S199 (memory problem)
#subjectList = ['S072', 'S075', 'S084', 'S117', 'S127', 'S128', 'S132', 'S133', 'S149', 'S183', 'S185', 
#               'S187', 'S191', 'S194', 'S195', 'S196', 'S197', 'S216'] 
subjectList = ['S218']
topchans = [31] #CHANGE AS NEEDED

wb = xlwt.Workbook()
sheet_ffr = wb.add_sheet('FFR_mag_plv_FCz')
style_bold = xlwt.easyxf('font: bold 1')

sheet_ffr.write(1, 0, 'Subject ID', style_bold)
sheet_ffr.write(0, 1, 'Mag_neural', style_bold) 
sheet_ffr.write(0, 3, 'Mag_not_purely_neural', style_bold) 
sheet_ffr.write(0, 5, 'Plv_neural', style_bold) 
sheet_ffr.write(0, 7, 'Plv_not_purely_neural', style_bold) 
sheet_ffr.write(0, 9, 'Mag_ABR', style_bold) 
sheet_ffr.write(1, 1, 'mag_masker_sum', style_bold) 
sheet_ffr.write(1, 2, 'mag_adpt_diff', style_bold)
sheet_ffr.write(1, 3, 'mag_masker_diff', style_bold)
sheet_ffr.write(1, 4, 'mag_prob_diff', style_bold)
sheet_ffr.write(1, 5, 'plv_masker_sum', style_bold) 
sheet_ffr.write(1, 6, 'plv_adpt_diff', style_bold)
sheet_ffr.write(1, 7, 'plv_masker_diff', style_bold)
sheet_ffr.write(1, 8, 'plv_prob_diff', style_bold)
sheet_ffr.write(1, 9, 'pos_abr', style_bold)
sheet_ffr.write(1, 10, 'neg_abr', style_bold)

wb.save(froot+'/FFR_mag_plv_FCz.xls')

fs = 16384
chop = 0.0e-3
chop2 = 50e-3

plvs = np.zeros(len(subjectList))
freqs = np.zeros(len(subjectList))
trialNums = np.zeros(len(subjectList))
fftMag = np.zeros(len(subjectList))
fftIndex = np.zeros(len(subjectList))

for k, subj in enumerate(subjectList):
    evoked_pos, evoked_neg, data_sumMasker, data_masker_diff, data_prob_diff, data_adpt_diff, numtrials_adpt = subjectProcessing(froot, subj, [1, 2], fs)
    diffWhole, diffProb, sumProb, diffMasked, sumMasked, diffMasker, sumMasker, diffAdpt, sumAdpt = diffSumEvd(evoked_pos, evoked_neg)
    PosMasker = evoked_pos[int(0.15*fs):int(0.25*fs)]
    NegMasker = evoked_neg[int(0.15*fs):int(0.25*fs)]
    peak_abr_pos, _ = fftPlot(PosMasker, fs, 0, 12e-3, 'ABR to positive masker')
    peak_abr_neg, _ = fftPlot(NegMasker, fs, 0, 12e-3, 'ABR to negative masker')
    peak_masker_sum, _ = fftPlot(sumMasker, fs, chop, chop2, 'MaskerSum')
    peak_adpt_diff, _ = fftPlot(diffAdpt, fs, chop, chop2, 'AdptDiff')
    peak_masker_diff, _ = fftPlot(diffMasker, fs, chop, chop2, 'MaskerDiff')
    peak_prob_diff, _ = fftPlot(diffProb, fs, chop, chop2, 'ProbDiff')
    
    plv_masker_sum, _, _ = plvAnalysis(data_sumMasker, 'MaskerSum')
    plv_adpt_diff, _, _ = plvAnalysis(data_adpt_diff, 'AdptDiff')
    plv_masker_diff, _, _ = plvAnalysis(data_masker_diff, 'MaskerDiff')
    plv_prob_diff, _, _ = plvAnalysis(data_prob_diff, 'ProbDiff')
    
    dictMat = {"mag_1k_sum": peak_masker_sum, "mag_500_adpt": peak_adpt_diff, "mag_500_75dB": peak_prob_diff, "mag_500_85dB": peak_masker_diff,
               "plv_masker_sum": plv_masker_sum, "plv_adpt_diff": plv_adpt_diff, "plv_masker_diff": plv_masker_diff, "plv_prob_diff": plv_prob_diff,
               "mag_abr_pos": peak_abr_pos, "mag_abr_neg": peak_abr_neg}
    savemat(froot+'/plv_mag/'+subj, dictMat)
    
    sheet_ffr.write(2+k, 0, subj)
    sheet_ffr.write(2+k, 1, str(peak_masker_sum)) 
    sheet_ffr.write(2+k, 2, str(peak_adpt_diff))
    sheet_ffr.write(2+k, 3, str(peak_masker_diff))
    sheet_ffr.write(2+k, 4, str(peak_prob_diff))
    sheet_ffr.write(2+k, 5, str(plv_masker_sum)) 
    sheet_ffr.write(2+k, 6, str(plv_adpt_diff))
    sheet_ffr.write(2+k, 7, str(plv_masker_diff))
    sheet_ffr.write(2+k, 8, str(plv_prob_diff))
    sheet_ffr.write(2+k, 9, str(peak_abr_pos))
    sheet_ffr.write(2+k, 10, str(peak_abr_neg))
    wb.save(froot+'/FFR_mag_plv_FCz.xls')
    
    #plvs[k] = plv_32_max
    #freqs[k] = f_max
    #trialNums[k] = numtrials
    #fftMag[k] = peak
    #index = index[0]
    #fftIndex[k] = index[0]
#dictDataArray = {"subjects": subjectList, "plvs": plvs, "freqs": freqs, "trialNums": trialNums, "fftMag": fftMag, "fftIndex": fftIndex}
#savemat('plvs_32_fftMags', dictDataArray)

#plv = spectral_connectivity(data_adpt_diff, method = 'plv')
#
##phi = np.asarray([phase475[10],  phase500[11], phase525[12]]) # check the index manually from the magnitude response
##f = np.asarray([475, 500, 525])
##delays = np.diff(np.unwrap(phi - 2*np.pi*f*chop)) * 1000./ (2*np.pi*np.diff(f))
##print delays
#
#dict500 = {"adpt500": diffAdpt, "evokedPos500": evoked_pos, "evokedNeg500": evoked_neg}
#savemat('evoekd500', dict500)