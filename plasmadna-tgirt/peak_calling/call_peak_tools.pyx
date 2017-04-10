import numpy as np


def merge_peaks(peak_start, peak_end):
    cdef:
        int tolerance_unprotected = 5
        int i = 0
        int j

    new_start = []
    new_end = []
    while i < len(peak_start)-2:
        new_start.append(peak_start[i])
        j = i
        while j < len(peak_start)-2 and peak_start[j+1] - peak_end[j] <= tolerance_unprotected:
            j += 1
        new_end.append(peak_end[j])
        j += 1
        i = j
    new_start.append(peak_start[i])
    new_end.append(peak_end[i])
    return np.array(new_start), np.array(new_end)

def maxSubArray(ls):
    '''
    #https://gist.github.com/alabid/3734606
    '''
    if len(ls) == 0:
       raise Exception("Array empty") # should be non-empty

    cdef:
        int runSum = ls[0]
        int maxSum = ls[0]
        int i = 0
        int start = 0
        int finish = 0
        int j

    for j in range(1, len(ls)):
        if ls[j] > (runSum + ls[j]):
            runSum = ls[j]
            i = j
        else:
            runSum += ls[j]

        if runSum > maxSum:
            maxSum = runSum
            start = i
            finish = j

    return start, finish
