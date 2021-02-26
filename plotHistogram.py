# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 12:51:12 2021

@author: oscbr226
"""

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from math import sqrt
from scipy.integrate import trapz 

def plotHistogram(histData, legends):
    '''Plots histograms and prints out statistics for multiple data series. 
    
    Plots histograms from multiple data series stored in histData, which is a 
    numpy array. The histogram is uses trapezoid based integration to give 
    a histogram to more resembles a line-based plot as compared to a normal
    histogram. The statistics that are printed out are the arithmetic mean,
    standard deviation, standard error of the mean, coefficient of variation
    and the number of data points in each dataset. Each data set also has an
    assoiciated legend to it. The function returns a figure and an axes object
    for further changes.
    
    Parameters 
    -----------
    histData : numpy object array with other numpy arrays with ints or floats
               Data to be plotted as histograms. Each index corresponds to 
               a dataset.
    legends : list of strings
              Strings associated with each dataset, providing information 
              where the data comes from.
              
    Returns
    -----------
    fig : figure object
          Contains figure that is created during the function. Is returned 
          so that the figure can be customized further.
    ax : axes objext
         Contains fig's axes object, which is created during the function.
         Returned for further customization options.
    '''
    
    fig, ax = plt.subplots(1,1)
    for i in range(len(histData)):
        data = histData[i]
        data = data[~np.isnan(data)]
        nrPoints = len(data)
        meandata = np.mean(data)
        stddata = np.std(data)
        semdata = stddata/sqrt(nrPoints)
        cvdata = stddata/meandata
        print(legends[i] + ':\nMean: '+ str(meandata) + '\nStd: ' +
              str(stddata) + '\nStd of mean: ' + str(semdata)
              + '\nCV: ' + str(cvdata) + '\nNr of data points: ' +
              str(nrPoints) + '\n')
        nrBins = np.round(sqrt(nrPoints))
        nrBins = nrBins.astype(int)
        [hist, binEdges] = np.histogram(data, nrBins)
        #Adding extra value to the histogram to get results similar to Matlab's
        #hist function.
        hist = np.append(hist, 0)
        ax.plot(binEdges, -hist/trapz(binEdges,hist), label=legends[i])
        
    ax.legend(loc=1)
    return fig, ax 
        
        