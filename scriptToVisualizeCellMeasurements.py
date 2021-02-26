# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 18:43:47 2021

@author: oscbr226
"""

import os
import numpy as np
from visualizeCellMeasurements import visualizeCellMeasurements

'''
Script to run visualizeCellMeasurements.py. This script is made so that it is
possible to visualize data generated from the Elf lab image analysis pipeline, 
which is run in Matlab. 

For this partiuclar application the cell cycle of E. coli is being studied to
answer questions regarding replication initiation control, how initiation 
is coupled to termination of replication and how this is connected to birth 
and division as well. Here we have performed time-lapse microscopy in micro-
fluidic devices where we have imaged a protein fused with the fluorescent 
protein venus. This fusion trails the replisome which allows us to track
the progression of replication faithfully. Furthermore, in the same strain
have have also labelled another protein with the fluorescent protein mCherry.
We have then inserted a sequence in the terminus region to which our second
fusion protein binds to which allows us to track the location of terminus. 
With this we can see where these components are located during the cell cycle
and draw conlusions whether they are colocalized or not.
'''

dataPath = os.getcwd()
lengthCutOff = 20
parentLengthCutOff = 10
daughtersLengthCutOff = 10
pixelSize = 11/150
heatMapThresh = 0.99
forkPlotBinScale = 25
erfSigmaCutOff = 1
axesLimits = {'xlim': [-3, 3], 'ylim': [4,1]}
titles = ['E. coli']
posInds = [np.arange(0,12)]
fluoChans = ['SeqA-venus', 'ParB-mCherry']

visualizeCellMeasurements(fluoChans, lengthCutOff, parentLengthCutOff,
                          daughtersLengthCutOff, dataPath, posInds, 
                          forkPlotBinScale, heatMapThresh, pixelSize, 
                          erfSigmaCutOff, titles, axesLimits)