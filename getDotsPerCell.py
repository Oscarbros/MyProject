# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 10:30:16 2021

@author: oscbr226
"""
import numpy as np

def getDotsPerCell(areas, counts, longs, lengths, minLength, pixelSize):
    '''Acquires number of dots per cell in area based bins.
    
    Function that acquires the number of fluorescent dots in cells within a 
    given cell size bin, where cell size corresponds to area. 
    
    Parameters
    -----------
    areas : numpy array with floats
            Array with areas of cells in frames where fluorescent dots were
            detected.
    counts : numpy array with floats
             Array with the number of dots in cells.  
    longs : numpy array with floats
            Array with internal long axis coordinates of fluorescent dots in
            cells.
    lengths : numpy array with floats
              Array with cell lengths in frames where fluorescent dots were 
              detected.
    minLength : float
                Parameter for the lower boundry of allowed internal long axis
                coordinates scaled with cell length.
    pixelSize : float
                Conversion factor between pixels and µm.
                
    Returns
    -----------
    dotsPerCell : numpy array with floats 
                  Normalized number of dots per cell within a given range of
                  cell sizes given in area.
    areaBins : numpy array with floats   
               Binned cell areas. 
    
    '''
    binScale = 50
    
    if not 'pixelSize' in locals():
        pixelSize = []
        
    if not 'minLength' in locals() or len(minLength) == 0:
        minLength = 0.44
        if pixelSize:
            minLength = minLength/pixelSize
            
    if pixelSize:
        binScale = binScale * pixelSize^2
        
    areaMinMax = np.quantile(areas, [0.005, 0.98])
    areaMin = np.min(areaMinMax)
    areaMax = np.max(areaMinMax)
    areaNrBins = np.rint((areaMax-areaMin)*binScale)
    areaNrBins = areaNrBins.astype(int)
    areaBins = np.linspace(areaMin, areaMax, areaNrBins)
    longs = longs - 0.5
    absLongLengths = np.abs(longs*lengths)
    
    dotsPerCell = np.zeros((areaNrBins-1,))
    for i in range(areaNrBins-1):
        selDots = np.logical_and(areas > areaBins[i], areas <= areaBins[i+1])
        selLongs = absLongLengths > minLength
        selCounts = counts[selDots]
        normFactor = np.sum(1/selCounts)
        dotsPerCell[i] = np.sum(selDots & selLongs)/normFactor
        
    areaBins = 0.5*(areaBins[0:-1]+areaBins[1:])
    
    return dotsPerCell, areaBins