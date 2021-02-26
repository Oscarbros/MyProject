# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 12:04:06 2021

@author: oscbr226
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plotDotLocations(longs, areas, lengths, counts, binScale, heatMapThresh,
                     meanBirthArea, meanDivisionArea, initiationArea):
    '''
    Plots a heatmap of fluorescent dots over the cell cycle.
    
    A function that creates a heatmap of fluorescent dots, which are sorted 
    based on size. In relation to this checkpoints in the cell cycle such as
    birth, division and initiation of replication can also be visualized. 
    

    Parameters
    ----------
    longs : numpy array with floats
            Array with internal long axis coordinates of fluorescent 
            dots in cells.
    areas : numpy array with floats
            Array with areas of cells in frames where fluorescent 
            dots were detected.
    lengths : numpy array with floats
              Array with cell lengths in frames where fluorescent dots 
              were detected.
    counts : numpy array with floats
             Array with the number of dots in cells.  
    binScale : int
               Scaling factor for bins in the heatmap. A larger number results 
               in more bins.
    heatMapThresh : float
                    Threshold for color values in the heatmap. Only allowed to 
                    be equal or larger than 0 and smaller than 1.
    meanBirthArea : float
                    Average birth area of the given cells. Only applicatple if 
                    the area has been calculated.
    meanDivisionArea : float
                       Average division area of the given cells. Only 
                       appplicable if the area has been calculated.
    initiationArea : float 
                     Average area at which cells initiate replication. Only 
                     applicable if the initiation area has been calculated and
                     replication has been tracked.

    Returns
    -------
    fig : figure object
          Figure returned so that it can be customized further.
    ax : axes object
         Axes object of the figure. Returned so that it can be modified 
         further.

    '''
    if len(longs) == 0 | len(areas) == 0 | len(lengths) == 0 | all(np.isnan(longs)):
        return
    
    if not 'binScale' in locals() or binScale == []:
        binScale = 20
    
    if not 'heatMapThresh' in locals() or heatMapThresh == []:
        heatMapThresh = 1
        
    if heatMapThresh <= 0 or heatMapThresh > 1:
        raise Exception('Parameter heatMapThresh must be in (0,1].')
        
    longs = longs - 0.5
    
    lmax = np.nanmax(lengths*longs)
    lmin = np.nanmin(lengths*longs)
    lMinMax = np.max([lmax, -lmin])
    lNrBins = np.rint(2*lMinMax*binScale)
    lNrBins = lNrBins.astype(int)
    if not lNrBins%2:
        lNrBins = lNrBins + 1
    lBins = np.linspace(-lMinMax, lMinMax, lNrBins)
    
    areaMinMax = np.quantile(areas, [0.005, 0.98])  
    areaMin = np.min(areaMinMax)
    areaMax = np.max(areaMinMax)
    areaNrBins = lNrBins
    areaBins = np.linspace(areaMin, areaMax, areaNrBins)
    
    heatMap = np.zeros((areaNrBins - 1, lNrBins-1))
    meanCellLengths = np.zeros((areaNrBins-1,))
    for i in range(areaNrBins-1):
        selDots = np.logical_and(areas > areaBins[i], areas <= areaBins[i+1])
        selLengths = lengths[selDots]
        selCounts = counts[selDots]
        selLong = longs[selDots]
        hist1, binEdges = np.histogram(selLong*selLengths, bins=lBins)
        normFactor = np.sum(1/selCounts)
        heatMap[i,:] = hist1/normFactor
        meanCellLengths[i] = np.mean(selLengths)
        
    if heatMapThresh < 1:
        heatMap[heatMap==0] = np.nan
        if meanBirthArea and meanDivisionArea:
            ind = np.logical_and(areaBins[0:-1] >= meanBirthArea, areaBins[1:] <= meanDivisionArea)
            data = heatMap[ind,:]
        else:
            data = heatMap
        threshold = np.quantile(data[:],heatMapThresh)
        heatMap[heatMap>threshold] = threshold
    
    fig, ax = plt.subplots(1,1)
    heatMap[np.isnan(heatMap)] = 0
    areaPlotRanges = 0.5*(areaBins[0:-1]+areaBins[1:])
    ax.imshow(heatMap, extent = [lBins[0], lBins[-1], areaPlotRanges[-1], areaPlotRanges[0]],  cmap=matplotlib.cm.jet)
    ax.plot(-0.5*meanCellLengths, 0.5*(areaBins[0:-1]+areaBins[1:]), color='w', lw=2)
    ax.plot(0.5*meanCellLengths, 0.5*(areaBins[0:-1]+areaBins[1:]), color='w', lw=2)
    ax.set_xlabel('Internal long axis coordinates')
    ax.set_ylabel(r'Area ($\mu$m$^2$)')
    
    if 'meanBirthArea' in locals() and meanBirthArea:
        ax.plot([lBins[0], lBins[-1]], [meanBirthArea, meanBirthArea], 'w--', lw=1.5)
    
    if 'meanDivisionArea' in locals() and meanDivisionArea:
        ax.plot([lBins[0], lBins[-1]], [meanDivisionArea, meanDivisionArea], 'w--', lw=1.5)
        
    if 'initiationArea' in locals() and initiationArea:
        ax.plot([lBins[0], lBins[-1]], [initiationArea, initiationArea], 'r--', lw=1.5)
        
    return fig, ax