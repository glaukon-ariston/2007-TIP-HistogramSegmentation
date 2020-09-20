#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Glaukon Ariston
# Date: 11.09.2020
# Abstract:
#       Matlab to Python transcription of the Matlab code from https://github.com/judelo/2007-TIP-HistogramSegmentation
#       Original Matlab code Copyright (c) 2016 Julie Delon
#
#       A non parametric approach for histogram segmentation_2006.pdf
#       https://scholar.google.fr/citations?user=X6MBHUUAAAAJ&hl=en
#       In this work, we propose a method to segment a 1-D histogram without a priori assumptions about the underlying density function.
#       https://github.com/judelo/2007-TIP-HistogramSegmentation
#       The following functions implement the Fine to Coarse Histogram Segmentation described in
#       J. Delon, A. Desolneux, J-L. Lisani and A-B. Petro, A non parametric approach for histogram segmentation
#       IEEE Transactions on Image Processing, vol.16, no 1, pp.253-261, Jan. 2007.
#       Usage :
#         u = double(imread('../images/lena.png'));
#         H = hist(u(:),0:255);
#         idx=FTC_Seg(H,0);
#       idx should contain the list of all minima separating the modes of H
#
# NOTES on Matlab to Python conversion
# https://realpython.com/matlab-vs-python/
# - for loops:
#       Matlab: for i=start:step:stop             -> stop is inclusive.
#       Python: for i in range(start,stop,step)   -> stop is exclusive.
# - Length of an array
#       Matlab: length(a)
#       Python: len(a), numpy.ndarray.size
# - Index into an array
#       Matlab: 1..length(a)
#       Python: 0..len(a)-1
#
# GNU Octave
#   > cd C:/app/dev/hp/projects/oss.uredjenazemlja.hr/kaptcha/papers/Threshold Selection/2007-TIP-HistogramSegmentation-master/code
# Modify code for computing local minima in FTC_Seg.m using findpeaks() as folows: 
#   % https://www.mathworks.com/matlabcentral/answers/44227-finding-local-minimums-maximums-for-a-set-of-data
#   Hinv = 1.01*max(H) - H;
#   [p,idx_min] = findpeaks(Hinv);
#
#   > pkg load signal
#   > u = double(imread('../images/lena.png'));
#   > H = hist(u(:),0:255);
#   > idx=FTC_Seg(H,0);
#
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import math

iterationNo = 0
def plot_histogram_modes(histogram, indices):
    global iterationNo
    valMax = np.amax(histogram)
    plt.plot(histogram)
    H = histogram
    N = histogram.size
    for i in range(indices.size):
        plt.axvline(x=indices[i], color='cyan', linestyle='dashed', alpha=0.5)
        plt.text(x=indices[i], y=valMax/2, s=str(i), alpha=0.7, color='#334f8d', rotation=90, va='center')
    plt.title("Histogram")
    plt.xlabel("Pixel Value")
    plt.ylabel("Pixel Count")
    plt.savefig('debug/%03d_histogram_modes_%d.png' % (iterationNo, indices.size))
    plt.close()
    iterationNo = iterationNo + 1


def plot_two_histograms(h1, h2):
    print(repr(h1))
    plt.clf()
    plt.style.use('seaborn-deep')
    bins = np.linspace(0, max(h1.size, h2.size), max(h1.size, h2.size))
    plt.hist(h1, bins, label='Interval', alpha=0.5)
    plt.hist(h2, bins, label='Isotonic Regression', alpha=0.5)
    plt.ylabel('Pixel count')
    plt.xlabel('Index')
    plt.title('Histogram Interval and its Isotonic Regression')
    plt.legend(loc='upper right')
    plt.show()


def plot_two_histograms1(h1, h2):
    print(repr(h1))
    plt.clf()
    plt.style.use('seaborn-deep')
    fig, ax = plt.subplot()

    x = np.arange(max(h1.size, h2.size))  # the label locations
    width = 0.1  # the width of the bars

    ax.bar(x - width/2, h1, width, label='Interval', alpha=0.5)
    ax.bar(x + width/2, h2, width, label='Isotonic Regression', alpha=0.5)

    ax.set_ylabel('Pixel count')
    ax.set_xlabel('Index')
    ax.set_title('Histogram Interval and its Isotonic Regression')
    ax.legend()
    fig.tight_layout()
    plt.show()


def pool_adjacent_violators(interval, nonDecreasing):
    '''
    https://en.wikipedia.org/wiki/Isotonic_regression
    Compute the isotonic regression of an histogram interval::np.array.
    nonDecreasing = boolean value saying if we want the non-decreasing (True) or decreasing regression (False)
    '''
    N = interval.size
    g = interval.copy().astype(np.float)
    if nonDecreasing:
        # non-decreasing regression
        for i in range(N-2, -1, -1):
            localSum = g[i]
            for j in range(i+1, N):
                if j == N-1 or g[j]*(j-i) >= localSum:
                    localSum = localSum / (j-i)
                    for k in range(i, j):
                        g[k] = localSum
                    break
                localSum = localSum + g[j]
    else:
        # decreasing regression
        for i in range(1, N):
            localSum = g[i]
            for j in range(i-1, -1, -1):
                if j == 0 or g[j]*(i-j) >= localSum:
                   localSum = localSum / (i-j);
                   for k in range(j+1, i+1):
                        g[k] = localSum
                   break
                localSum = localSum + g[j]
    return g


def entropy(x, y):
    '''
    Function computing the entropy between x and y.
    x and y must be in the interval [0,1]
    '''
    v = -math.log10(1-y) if x == 0.0 else (
            -math.log10(y) if x == 1.0
            else (x*math.log10(x/y) + (1.0-x)*math.log10((1.0-x)/(1.0-y)))
        )
    return v



def max_entropy(h, a, b, e, nonDecreasing):
    '''
    Compute the maximum entropy of the histogram h(a:b) for the increasing or decreasing hypothesis
    nonDecreasing = boolean value indicating if we test the increasing or decreasing hypothesis
    h = histogram
    e = parameter used to compute the entropy
    '''
    g = h[a:b+1]
    L = g.size
    decreas = pool_adjacent_violators(g, nonDecreasing)
    #plot_two_histograms(g, decreas)

    # integrate signals
    g = np.cumsum(g)
    decreas = np.cumsum(decreas)

    # meaningfullness threshold
    N = g[L-1]
    threshold = (math.log(L*(L+1)/2) + e*math.log(10)) / N

    # search the most meaningfull segment (gap or mode)
    maxEntropy = 0.0
    for i in range(0, L):
        for j in range(i, L):
            r = g[j] if i == 0 else g[j] - g[i-1]
            r = r / N
            p = decreas[j] if i == 0 else decreas[j] - decreas[i-1]
            p = p / N
            v = entropy(r, p)
            if v > maxEntropy:
                maxEntropy = v
    maxEntropy = (maxEntropy - threshold) * N
    return maxEntropy


def fine_to_coarse_histogram_segmentation(H, e):
    '''
    Implements the Fine to Coarse Histogram Segmentation described in
    J. Delon, A. Desolneux, J-L. Lisani and A-B. Petro, A non parametric approach for histogram segmentation
    IEEE Transactions on Image Processing, vol.16, no 1, pp.253-261, Jan. 2007.
    Usage:
        u = double(imread('../images/lena.png'));
        H = hist(u(:),0:255);
        idx=FTC_Seg(H,0);
    idx should contain the list of all minima separating the modes of H

    FTC_seg
    H = the integer-valued histogram to be segmented.
    e = parameter of the segmentation (corresponds to e = -log10(epsilon) in the paper)
    large e => coarse segmentation
    small e => fine segmentation
    '''
    N = H.size

    # find the list of local minima and maxima of H
    idxs_max, _ = signal.find_peaks(H)
    idxs_min, _ = signal.find_peaks(-H)
    idxs = np.sort(np.concatenate([np.array([0]), idxs_min, idxs_max, np.array([N-1])]))

    # find if idxs starts with a minimum or a maximum
    begins_with_min = H[idxs[0]] < H[idxs[1]]

    # FILL THE LIST OF ENTROPIES FOR ALL MODES
    # The merging of two contiguous modes [a,b] and [b,c] can be done in two ways,
    # either by using the maximum M1 on [a,b] and by testing the decreasing hypothesis on [M1,c],
    # or by using the maximum M2 on [b,c] and by testing the increasing hypothesis on [a,M2].
    # For each configuration, we compute the entropy of the worst interval against the considered hypothesis.
    K = idxs.size
    maxEntropies = np.zeros((K-3,))
    # Loop on all optimas
    for k in range(0, K-3):
        # decide if we want to test the increasing or decreasing hypothesis on [idxs(k), idxs(k+3)]
        nonDecreasing = 1 if (not begins_with_min and k % 2 == 1) or (begins_with_min and k % 2 == 0) else 0
        # compute the max entropy on the interval [k,k+3]
        maxEntropies[k] = max_entropy(H, idxs[k], idxs[k+3], e, nonDecreasing)

    # MERGING of MODES
    # [idxs(kmin), idxs(kmin+3)] is the first interval to merge
    kmin = np.argmin(maxEntropies)
    valmin = maxEntropies[kmin]

    while maxEntropies.size > 0 and valmin < 0:
        plot_histogram_modes(H, idxs)
        # update the list of min, max
        idxs = np.concatenate([idxs[0:kmin+1], idxs[kmin+3:]])
        maxEntropies = np.concatenate([maxEntropies[0:kmin+1], maxEntropies[kmin+3:]])      # No need to worry about upper bound out of range in Python.
        maxEntropies = maxEntropies[0:idxs.size-3+1]

        # update max_entropy around the removed optima
        for j in range(max(kmin-2,0), min(kmin,maxEntropies.size)):
            # decide if increasing or decreasing
            nonDecreasing = 1 if (begins_with_min and j % 2 == 1) or (not begins_with_min and j % 2 == 0) else 0
            # update the max entropy on the interval [k,k+3]
            maxEntropies[j] = max_entropy(H, idxs[j], idxs[j+3], e, nonDecreasing)
        if maxEntropies.size > 0:
            kmin = np.argmin(maxEntropies)
            valmin = maxEntropies[kmin]

    idxs = idxs[0::2] if begins_with_min else idxs[1::2]
    return idxs
