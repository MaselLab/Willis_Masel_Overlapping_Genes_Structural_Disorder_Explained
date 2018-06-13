This is the main linear model used in this project. It is a linear mixed model with two fixed effects (Designation and Frame) and two random effects (Homology Group and Species). 

The model comes with a data file to be read in. To do this, the user needs to specify the path to the data file in the script. 

The data file contains the ancestral and novel overlapping genes for which relative ages could be determined. The ISD values associated with each overlapping gene are the ISD values for the overlapping section only, and not the entire gene. The controls included in this file are the artificially-frameshifted non-overlapping viral controls.

The script is annotated for ease of use. 