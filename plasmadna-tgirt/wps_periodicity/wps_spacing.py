s
#!/usr/bin/env python

import numpy as np
import pandas as pd
from statsmodels.tsa.filters.filtertools import recursive_filter
from scipy.signal import detrend, periodogram
from rpy2 import robjects
from rpy2.robjects import numpy2ri

periodogram = robjects.r('spec.pgram')
def r_spectrum(signal):
    filtered_signal = detrend(signal, type='linear')
    filtered_signal = recursive_filter_function(filtered_signal)
    filtered_signal = demean(filtered_signal)
    res = periodogram(numpy2ri.numpy2ri(filtered_signal),
                        pad=0.3,tap=0.3,
                        span=2,plot=False,detrend=True,
                        demean=True)
    freq = numpy2ri.ri2py(res.rx2('freq'))
    psd = numpy2ri.ri2py(res.rx2('spec'))
    return 1/freq, psd


def shift_array(signal):
    signal = np.array(signal)
    signal = np.append(signal[:300], signal)
    return signal


filter_freq = 1.0/np.arange(5,100,4)
def recursive_filter_function(signal):
    signal = shift_array(signal)
    filter_length = len(filter_freq)
    signal = recursive_filter(signal, filter_freq)[300:]
    return signal


def demean(arr):
    return arr-arr.mean()


def weighted_smooth(w):                        
        def smoothing(arr):
            return (arr*w).mean()
        return smoothing


def daniell_smoother(arr):
    arr = pd.Series(arr)\
        .rolling(window = 5)\
        .apply(weighted_smooth(np.array([1,1,1,1,0.5])))
    return np.array(arr)

def daniell_spectrum(signal):
    # detrending and smoothing before spectrogram
    filtered_signal = recursive_filter_function(signal)
    filtered_signal = detrend(filtered_signal, type='linear')
    filtered_signal = demean(filtered_signal)
    freq, psd = periodogram(filtered_signal, detrend='linear')
    psd = daniell_smoother(psd)
    return 1/freq, psd
