# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 10:56:07 2021

@author: oscbr226
"""

import numpy as np
from getDotsPerCell import getDotsPerCell
from scipy.optimize import curve_fit 
from scipy.special import erf
from math import sqrt
import matplotlib 
import matplotlib.pyplot as plt

def plotDotsPerCell(areas, counts, longs, lengths, minLength, avgDivisionArea,
                    pixelSize, showStatistics):
    '''Plots normalized fluorescent dots per cell in area bins.
    
    Function that plots normalized fluorescent dots per cell at different 
    area ranges. The areas correspond to different stages from birth to 
    division. Furthermore, the average initiation areas of all cells is 
    determined.
    
    Parameters
    ------------
    areas : numpy array with floats
            Array with areas of cells in frames where fluorescent 
            dots were detected.
    counts : numpy array with floats
             Array with the number of dots in cells.  
    longs : numpy array with floats 
            Array with internal long axis coordinates of fluorescent 
            dots in cells.
    lengths : numpy array with floats
              Array with cell lengths in frames where fluorescent dots 
              were detected.
    minLength : float (optional)
                Parameter for the lower boundry of allowed internal long 
                axis coordinates scaled with cell length.
    avgDivisionArea : float 
                      Average area at which cells divide.
        
    pixelSize : float (optional)
                Conversion factor between pixels and µm.
    showStatistics : Boolean
                     Flag to determine whether initiation area statistics
                     should be printed or not.
                      
    Returns 
    -----------
    dotsPerCell : numpy array with floats.
                  Normalized number of dots per cell within a given range 
                  of cell sizes given in area.
    areaBins : numpy array with floats   
               Binned cell areas.
    theFit : Tuple containing numpy arrays containing floats. 
             The first index in the tuple contains the constants 
             found by curve_fit when fitting dots per cell and areas 
             to errorfunc. The second index contains a 2D numpy array
             with estimated covariance of the constants in the first
             index.
    fig : figure object
          Contains figure that is created during the function. 
          Is returned so that the figure can be customized further.
    ax : axes objext
         Contains fig's axes object, which is created during the 
         function. Returned for further customization options.'''
    
    if np.isnan(avgDivisionArea):
        dotsPerCell = []
        areaBins = []
        theFit = []
        return dotsPerCell, areaBins, theFit
    
    if not 'pixelSize' in locals():
        pixelSize = []
        
    if not 'minLength' in locals():
        minLength = []
    
    if not 'showStatistics' in locals():
        showStatistics = []
        
    startPoint = np.array([2, 2, 0.2])
    if pixelSize:
        startPoint[1:] = startPoint[1:]/(pixelSize**2)
    
    [dotsPerCell, areaBins] = getDotsPerCell(areas, counts, longs, lengths, minLength, pixelSize)
    
    areaBins = areaBins[~np.isnan(dotsPerCell)]
    dotsPerCell = dotsPerCell[~np.isnan(dotsPerCell)]
    
    selAreaBins = areaBins < avgDivisionArea
    areaBinsFit = areaBins[selAreaBins]
    dotsPerCellFit = dotsPerCell[selAreaBins]
    theFit = curve_fit(errorfunc,areaBinsFit,dotsPerCellFit)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(areaBins, dotsPerCell, 'bs')
    ax.plot(areaBinsFit,errorfunc(areaBinsFit, theFit[0][0],theFit[0][1],theFit[0][2]), color='r')
    ylims = ax.get_ylim()
    ax.plot([theFit[0][1], theFit[0][1]],[0, ylims[1]], 'k--',  lw = 2)
    ax.set_ylim([0, ylims[1]])
    ax.set_xlabel(r'Area ($\mu$m$^2$)')
    ax.set_ylabel('Dots per cell')
    
    
    if showStatistics:
        print('*********Statistics****************')
        print('Initiation volume: ' + str(theFit[0][1]))
        print('Std initaition volume: ' + str(theFit[0][2]))
        print('***********************************')
       
    return dotsPerCell, areaBins, theFit, fig, ax
    
def errorfunc(x, b, x0, s):
    '''Function for fitting dots per cell and areas to.
    
    Function used for error-function-based-fitting to acquire the average
    initiation area for all cells in the dataset.
    
    Parameters 
    -----------
    x : array or float
        Array with cell areas
    b : float
        First degree constant for dots per cells.
    x0 : float  
         Constant corresponding to averge initiation area.
    s : float 
        Constant corresponding to spread of initiation.
    '''
    return b*(1+erf((x-x0)/(sqrt(2)*s)))
    