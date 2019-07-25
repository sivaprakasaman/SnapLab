#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:51:03 2019

@author: parida
"""
import numpy as np


def numel(inVal):
    if isinstance(inVal, (list, str)):
        outLen= len(inVal)
    elif isinstance(inVal, (int, int, float)):
        outLen= 1
    else: # numpy array and others?
        outLen= len(inVal)
    
    return outLen

def dbspl(inSignal):
    pRef= 20e-6
    outSignal= 20*np.log10(inSignal/pRef)
    
    return outSignal

def calc_dbspl(inSignal):
    pRef= 20e-6
    splVal= 20*np.log10(rms(inSignal)/pRef)
    
    return splVal

def rms(inSignal):
    rms_val = np.sqrt(np.mean(inSignal**2))
    return rms_val