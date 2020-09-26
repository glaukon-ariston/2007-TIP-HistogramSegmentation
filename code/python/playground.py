#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Glaukon Ariston
# Date: 26.09.2020
# Abstract:
#   Debugging harness for histogramSegmentation.fine_to_coarse_histogram_segmentation(H, e=0)
#
# (findpeaks in Python get different results than in MatLab)
# print('idx = [%s];' % (' '.join(['%d' % (v+1,) for v in idxs]),))
# idx = [1 48 49 50 67 68 69 71 72 73 74 76 77 79 81 82 83 84 85 92 93 100 102 103 110 112 113 114 116 119 120 127 128 130 136 137 138 145 147 148 150 155 169 170 171 173 176 178 186 188 189 196 197 199 200 202 203 206 207 209 211 212 221 222 234 235 237 239 241 243 244 245 256];
#
# print('val[%d] = [%s]' % (len(maxEntropies), ' '.join(['%f' % (v,) for v in maxEntropies])))
# val[70] = [-7.044258 -5.163717 469.061042 -5.331204 -2.602968 -1.984528 -2.600699 -2.296912 -2.674793 -1.536530 -2.815287 -2.116291 -2.991723 -2.703314 -1.531241 -0.981140 -3.991059 7.867301 -4.869646 69.468887 -2.245196 -4.134025 41.480205 -4.046769 -2.674656 -2.402819 -2.774352 0.178781 -4.328724 47.723254 -4.186474 -3.127704 19.860393 -3.645382 -3.957270 17.571696 -1.331898 -2.549457 -3.221580 12.974305 453.162900 -4.956052 -2.065192 -0.415594 -1.134470 -4.444906 107.802875 -3.603364 -3.993935 20.479040 -2.063956 -2.681444 -2.386676 -1.899041 -2.988679 -0.823756 -3.226691 -0.842533 -0.517108 -3.725309 223.395585 -5.520483 170.863611 -4.879836 0.966133 -3.032572 -2.957388 -2.141432 -2.179776 -4.337235];
# val[70] = [-7.044258 -5.163717 469.061042 -5.331204 -2.602968 -1.984528 -2.600699 -2.296912 -2.674793 -1.536530 -2.815287 -2.116291 -2.991723 -2.703314 -1.531241 -0.981140 -3.991059 7.867301 -4.869646 69.468887 -2.245196 -4.134025 41.480205 -4.046769 -2.674656 -2.402819 -2.774352 0.178781 -4.328724 47.723254 -4.186474 -3.127704 19.860393 -3.645382 -3.957270 17.571696 -1.331898 -2.549457 -3.221580 12.974305 453.162900 -4.956052 -2.065192 -0.415594 -1.134470 -4.444906 107.802875 -3.603364 -3.993935 20.479040 -2.063956 -2.681444 -2.386676 -1.899041 -2.988679 -0.823756 -3.226691 -0.842533 -0.517108 -3.725309 223.395585 -5.520483 170.863611 -4.879836 0.966133 -3.032572 -2.957388 -2.141432 -2.179776 -4.337235];
# val[68] = [-7.044258 -5.331204 -2.602968 -1.984528 -2.600699 -2.296912 -2.674793 -1.536530 -2.815287 -2.116291 -2.991723 -2.703314 -1.531241 -0.981140 -3.991059 7.867301 -4.869646 69.468887 -2.245196 -4.134025 41.480205 -4.046769 -2.674656 -2.402819 -2.774352 0.178781 -4.328724 47.723254 -4.186474 -3.127704 19.860393 -3.645382 -3.957270 17.571696 -1.331898 -2.549457 -3.221580 12.974305 453.162900 -4.956052 -2.065192 -0.415594 -1.134470 -4.444906 107.802875 -3.603364 -3.993935 20.479040 -2.063956 -2.681444 -2.386676 -1.899041 -2.988679 -0.823756 -3.226691 -0.842533 -0.517108 -3.725309 223.395585 -5.520483 170.863611 -4.879836 0.966133 -3.032572 -2.957388 -2.141432 -2.179776 -4.337235];
#
# GNU Octave
# logId = fopen('val_matlab.m', 'w');
# fprintf(logId, "val[%d] = [%s]\n", length(val), sprintf("%f ", val))
# fprintf(logId, 'j %d inc %d val(j) %f\n', j, inc, val(j))
# val[68] = [-7.044258 -5.331204 -2.602968 -1.984528 -2.600699 -2.296912 -2.674793 -1.536530 -2.815 287 -2.116291 -2.991723 -2.703314 -1.531241 -0.981140 -3.991059 7.867301 -4.869646 69.468887 -2.2 45196 -4.134025 41.480205 -4.046769 -2.674656 -2.402819 -2.774352 0.178781 -4.328724 47.723254 -4 .186474 -3.127704 19.860393 -3.645382 -3.957270 17.571696 -1.331898 -2.549457 -3.221580 12.974305 453.162900 -4.956052 -2.065192 -0.415594 -1.134470 -4.444906 107.802875 -3.603364 -3.993935 20.4 79040 -2.063956 -2.681444 -2.386676 -1.899041 -2.988679 -0.823756 -3.226691 -0.842533 -0.517108 - 3.725309 223.395585 -5.520483 170.863611 -4.879836 0.966133 -3.032572 -2.957388 -2.141432 -2.1797 76 -4.337235 ]

import imageio
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters
from scipy import signal
import histogramSegmentation
import math
import os.path
import os
import glob


def playgroundHistogramSegmentation(imagePath):
    img = np.array(Image.open(imagePath).convert('L'))
    #histogram, bins = np.histogram(img.ravel(), 256, [0,256])
    histogram  = np.bincount(img.ravel(), minlength=256)

    plt.plot(histogram)

    H = histogram
    N = histogram.size
    #idxs_max, _ = signal.find_peaks(H)
    #idxs_min, _ = signal.find_peaks(-H)
    #idxs = np.sort(np.concatenate([np.array([0]), idxs_min, idxs_max, np.array([N-1])]))
    #for i in range(idxs_max.size):
    #    plt.axvline(x=idxs_max[i], color='cyan')
    #for i in range(idxs_min.size):
    #    plt.axvline(x=idxs_min[i], color='magenta')

    indices = histogramSegmentation.fine_to_coarse_histogram_segmentation(H, e=0)
    for i in range(indices.size):
        plt.axvline(x=indices[i], color='cyan')

    #_ = plt.hist(img.ravel(), 256, [0,256])
    #plt.imshow(img, cmap='gray', interpolation='none')
    plt.title("Histogram")
    plt.show()


def main():
    playgroundHistogramSegmentation('../../images/lena.png')


if __name__ == "__main__":
  main()
