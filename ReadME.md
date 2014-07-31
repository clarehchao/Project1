Evaluation of Positron Range Spread for Positron Emission Tomography (PET) Imaging
===================================================================================================

Purpose
-------------------
The python script was written to process the positron position data tracked using GEANT4 Monte Carlo simulations.
The goal is to determine the effect of magnetic field in positron range of long-range positron emitters (Ga-68, I-124, or Rb-82) in positron emission tomography (PET) imaging.
Monte Carlo simulations were developed to track positron in some medium of interaction, e.g. water, tissue, or air.
The positron range spread is evaluated based on the position where positron annhilation occurs.


Details
--------------------
- Read in the txt files from the simulations
- Discretize the positron position in a volume of a known voxel dimension
- Evaluate the line profile of a image slice of the constructed volume
- Applied curve fitting to the line profile to determine full width half maximum (FWHM) and full width tenth maximum (FWTM)
- Utilized Numpy and Scipy packages in Python (version 2.7.3)

Example Results
--------------------
The result illustrates the effect of 7-Tesla magnetic field on the positron range spread in two image directions (in-plane and transverse)

### In-plane direction
-----------------------
This figure shows a in-plane slice image of the positron position volume 
![fig1](https://github.com/clarehchao/Project1/blob/master/data/XPositronPostPosition_Run1_50.jpg "A slice image of the positron position volume")
The below figure is the line profile of the in-plane image slice of the positron position volume.
![fig2](https://github.com/clarehchao/Project1/blob/master/data/XLineProfileFit_PostPosition_Run1_50.jpg "Line profile in X-diretion of the positron position volume")

### Transverse direction
-----------------------
The figure shows a transverse slice image of the positron position volume.
![fig3](https://github.com/clarehchao/Project1/blob/master/data/ZPositronPostPosition_Run1_50.jpg "A slice image of the positron position volume")
The below figure is the line profile and of the transverse image slice of the positron position volume.
![fig4](https://github.com/clarehchao/Project1/blob/master/data/ZLineProfileFit_PostPosition_Run1_50.jpg "Line profile in X-diretion of the positron position volume")













