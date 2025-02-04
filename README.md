# MyProject
Repository for my project in the course Advanced Scientific Programming.
In the project I will rework some Matlab code where I create plots with data from .mat files. 
The plots include heat maps, fitting an error function to data, distributions of data and 
plots with errorbars. The plotting is done on data generated by post-processing calculations
from our image analysis pipeline.

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

The file to run all the functions is called scriptToRunVisualizeCellMeasurements.py.
