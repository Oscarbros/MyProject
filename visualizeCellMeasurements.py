# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 15:36:36 2021

@author: oscbr226
"""

import numpy as np 
from scipy.io import loadmat
import os
from math import sqrt
from math import log
from plotDotsPerCell import plotDotsPerCell
from plotDotLocations import plotDotLocations
from plotHistogram import plotHistogram
import matplotlib
import matplotlib.pyplot as plt

def visualizeCellMeasurements(fluoChans, lengthCutOff, parentLengthCutOff, daughtersLengthCutOff,
                              dataPath, strainPosInds, binScale, heatMapThresh, pixelSize,
                              erfSigmaCutOff, titles, axesLimits):
    '''
    Function that visualizes different cell measurements.
    
    Used for visualization of image analysis output from the Elf lab 
    image analysis pipeline and post pipline processing. This data is stored in 
    a .mat file. The function visualizes the numbers of dots in cells sorted on 
    size as well as the distribution of the dots over the cell cycle. 
    Furthermore, important values such as birth area, division area, growth 
    rate generation time and initiation area are also calculated and 
    visualized. These can provide insight how strains differ between one 
    another when it comes to the cell cycle.

    Parameters
    ----------
    fluoChans : list of strings
                Contains information regarding the fluorophore used in the 
                experiment.
    lengthCutOff : int
                   Cut-off value for how many frames a cell must have been 
                   detected to be used in the analysis. Here it is used for 
                   loading data.
    parentLengthCutOff : int
                         Cut-off value for how many frames a cell's parent
                         must have been detected to be used in the analysis.
                         Here it is used for loading data.
    daughtersLengthCutOff : int
                            Cut-off value for how many frames a cell's 
                            daughters must have been detected to be used in
                            the analysis. Here it is used for loading data.
    dataPath : string
               String containing the path from where to load data from a .mat
               file.
    strainPosInds : numpy array with ints
                    Numpy array with indices used to differentate between
                    different strains. Its primary use is when data from 
                    multiple strains is store in the same .mat file.
    binScale : int
               Scaling factor for bins in the heatmap. A larger number results 
               in more bins.
    heatMapThresh : float
                    Threshold for color values in the heatmap. Only allowed to 
                    be equal or larger than 0 and smaller than 1.
    pixelSize : float
                Conversion factor between pixels and um.
    erfSigmaCutOff : float
                     Cut-off value for the standard deviation of the average
                     initiation area. If the standard deviation exceeds or is
                     equal to this value the fit is considered insignificant, 
                     meaning that no initiation area could be determined.
    titles : list of strings
             Contains information about strains, genotypes or other useful 
             information to highlight what is plotted.
    axesLimits : dictonary with lists
                 Dictionary that contains limits on the x- and y-axis of the 
                 heatmaps.

    Returns
    -------
    None.

    '''
    nrFluoChans = len(fluoChans)
    nrStrains = len(strainPosInds)
    
    matfileData = loadmat(os.path.join(dataPath, 'measuredCells_' + str(lengthCutOff) + 
                        '_' + str(parentLengthCutOff) + '_' + str(daughtersLengthCutOff) 
                        + '.mat'))
    
    allPosLongs = matfileData.get('allPosLongs')
    allPosAreas = matfileData.get('allPosAreas')
    allPosLengths = matfileData.get('allPosLengths')
    allPosCounts = matfileData.get('allPosCounts')
    allBirthAreas = matfileData.get('allBirthAreas')
    allDivisionAreas = matfileData.get('allDivisionAreas')
    allGrowthRates = matfileData.get('allGrowthRates')
    allGenerationTimes = matfileData.get('allGenerationTimes')
    
    xlims = axesLimits.get('xlim')
    ylims = axesLimits.get('ylim')
        
    grForPlots = []
    grLegends = []
    gtForPlots = []
    gtLegends = []
    ind = 0
    
    for i in range(nrStrains):
        
        posInds = strainPosInds[i]
        bAreas = allBirthAreas[posInds]
        bAreas = unpackRawCellStats(bAreas)*pixelSize**2
        dAreas = allDivisionAreas[posInds]
        dAreas = unpackRawCellStats(dAreas)*pixelSize**2
        gRates = allGrowthRates[posInds]
        gRates = unpackRawCellStats(gRates)
        gTimes = allGenerationTimes[posInds]
        gTimes = unpackRawCellStats(gTimes)
        
        nrBirths = len(bAreas)
        meanBArea = np.nanmean(bAreas)
        stdBArea = np.nanstd(bAreas)
        semBArea = stdBArea/sqrt(nrBirths)
        cvBArea = stdBArea/meanBArea
        
        nrDivisions = len(dAreas)
        meanDArea = np.nanmean(dAreas)
        stdDArea = np.nanstd(dAreas)
        semDArea = stdDArea/sqrt(nrDivisions)
        cvDArea = stdDArea/meanDArea
        
        for j in range(nrFluoChans):
            
            longs = allPosLongs[posInds] 
            longs = unpackDotCellStats(longs,j)
            areas = allPosAreas[posInds]
            areas = unpackDotCellStats(areas,j)*pixelSize**2
            lengths = allPosLengths[posInds]
            lengths = unpackDotCellStats(lengths,j)*pixelSize
            counts = allPosCounts[posInds]
            counts = unpackDotCellStats(counts,j)
            
            _, _, theFit, currFig, currAx = plotDotsPerCell(areas, counts, longs, lengths, [], meanDArea, [], [])
            currAx.set_title(titles[i] + ' ' + fluoChans[j])
            
            if np.logical_and(len(theFit) != 0,theFit[0][2] < erfSigmaCutOff):
                meanInitArea = theFit[0][1]
            else:
                meanInitArea = []
                currAx.remove()
                
            currFig, currAx = plotDotLocations(longs, areas, lengths, counts, binScale, heatMapThresh, meanBArea, meanDArea, meanInitArea)
            currAx.set_xlim(xlims)
            currAx.set_ylim(ylims)
            currAx.set_title(titles[i] + ' ' + fluoChans[j])
            
            if meanInitArea:
                stdIArea = theFit[0][2]
                nonNanInds = ~np.isnan(longs)
                nonNanLongs = longs[nonNanInds]
                nrIAreas = len(nonNanLongs)
                semIArea = stdIArea/sqrt(nrIAreas)
                cvIArea = stdIArea/meanInitArea
                print(titles[i] + ':\nMean BA: ' + str(meanBArea) + ', Std BA: ' +
                      str(stdBArea) + ', SEM BA: ' +  str(semBArea) + 
                      ', CV BA: ' + str(cvBArea) + ', Nr births: ' 
                      + str(nrBirths) + '\nMean DA: ' + str(meanDArea) + 
                      ', Std DA: ' + str(stdDArea) + ', SEM DA: ' 
                      + str(semDArea) + ', CV DA: ' + str(cvDArea) 
                      + ', Nr divisions: ' + str(nrDivisions) + 'n\Mean IA: ' 
                      + str(meanInitArea) + ', Std IA: ' + str(stdIArea) + 
                      ', SEM IA: ' + str(semIArea) + ', CV IA: ' 
                      + str(cvIArea) + ', Nr dot locations: ' 
                      + str(nrIAreas))
            else:
                print(titles[i] + ':\nMean BA: ' + str(meanBArea) + ', Std BA: ' +
                      str(stdBArea) + ', SEM BA: ' +  str(semBArea) + 
                      ', CV BA: ' + str(cvBArea) + ', Nr births: ' 
                      + str(nrBirths) + '\nMean DA: ' + str(meanDArea) + 
                      ', Std DA: ' + str(stdDArea) + ', SEM DA: ' 
                      + str(semDArea) + ', CV DA: ' + str(cvDArea) 
                      + ', Nr divisions: ' + str(nrDivisions))
                   
        growthRateGenTimeBased = 60*log(2)/gTimes
        genTimesGrowthRateBased = log(2)/gRates        
        
        gRates = gRates*60
        
        bothFluorTypes = ' '.join(fluoChans)
        grForPlots.append(gRates)
        grLegends.append(titles[i] + ' ' + bothFluorTypes + '\nExponential fit based')
        grForPlots.append(growthRateGenTimeBased) 
        grLegends.append(titles[i] + ' ' + bothFluorTypes + '\nDetection time based')
        ind += 1
        gtForPlots.append(gTimes)
        gtLegends.append(titles[i] + ' ' + bothFluorTypes + '\nDetection time based')
        gtForPlots.append(genTimesGrowthRateBased) 
        gtLegends.append(titles[i] + ' ' + bothFluorTypes + '\nExponential fit based')
        
    print('*******Growth rates*******')
    currFig, currAx = plotHistogram(grForPlots, grLegends)
    currAx.set_xlim([0, 2.5])
    ylims = currAx.get_ylim()
    currAx.set_ylim([0,ylims[1]])
    currAx.set_xlabel(r'Growth rate (h$^{-1}$)')
    
    print('*******Generation times*******')
    currFig, currAx = plotHistogram(gtForPlots, gtLegends)
    currAx.set_xlim([15, 120])
    ylims = currAx.get_ylim()
    currAx.set_ylim([0,ylims[1]])
    currAx.set_xlabel('Generation time (min)')
    
    plt.show()    
    
def unpackDotCellStats(cellMeasArray,fluoChanInd):
    '''
    Function to unpack specific storage format from .mat file.
    
    Function that unpacks data from a .mat file that is stored into a cell 
    array where each cell contains another cell array where each cell contains
    an array with doubles. Each cell in the first layer corresponds to one 
    image position while the inner layer corresponds to one fluorescence 
    channel.
    
    Parameters
    ----------
    cellMeasArray : numpy array
                    Numpy array that is structured according to the function
                    loadmat in the scipy.io package.
    fluoChanInd : int
                  Index for the corresponding fluorescence channel.

    Returns
    -------
    unpackedStats : numpy array
                    Numpy array with unpacked and concatenated cells stats.

    '''
    nShape = np.shape(cellMeasArray)
    unpackedStats = []
    for i in range(nShape[0]):
        unpackedStats = np.concatenate((unpackedStats,cellMeasArray[i][0][fluoChanInd][0][0]))
    return unpackedStats
                 

def unpackRawCellStats(cellStats):
    '''
    Function to unpack specific storage format from .mat file.
    
    Function that unpacks data from a .mat file that is stored into a cell 
    array where each cell contains an array with doubles and then concatenates
    each array into a numpy array. Each array with doubles corresponds to data
    from one position.
    
    Parameters
    ----------
    cellMeasArray : numpy array
                    Numpy array that is structured according to the function
                    loadmat in the scipy.io package.
                    
    Returns
    -------
    unpackedStats : numpy array
                    Numpy array with unpacked and concatenated cells stats.

    '''
    nShape = np.shape(cellStats)
    unpackedStats = []
    for i in range(nShape[0]):
        unpackedStats = np.concatenate((unpackedStats, cellStats[i][0][0]))
        
    return unpackedStats
        
        
        
        
        
        