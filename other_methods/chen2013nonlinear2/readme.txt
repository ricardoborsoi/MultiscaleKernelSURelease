This demo box is related to nonlinear unmixing algorithm KHype and its spatial regularized version


Version 11-Mar.-2013
If you have any questions, please contact:
chen@uncie.fr
or cedric.richard@unice.fr

----- Executable files -----
Demo_KHype:       The nonlinear unmixing algorithm KHype in the following paper:

J. CHEN, C. RICHARD and P. HONEINE, 
"Nonlinear unmixing of hyperspectral data based on a linear-mixture/nonlinear-fluctuation model". 
IEEE Transactions on signal processing, 2013.

Demo_Spatial_KHype:  Spatial regularized nonlinear unmixing in the following paper:

J.CHEN, C. RICHARD  and P. HONEINE, 
"Nonlinear estimation of material abundances of hyperspectral images with l1-norm spatial regularization". 
IEEE Trans. on Geoscience. (RQ).

----- Complementary fucntions -----
(For KHype:)
AbundanceGen:     Abundance generation function;
generate_image:   Test image generation function;
hypermix:         Mixing endmembers with models;
ErrComput:        Computation of Estimation error and standard deviation;

(For Spatial regualarized KHype:)
grid_image:       Generate IM1 in the paper
ConvH, ConvHT, IpHHT:  functions to calculate H*(.) , (.)*H', (I+HH')

----- Database -------
endmembers.mat:   The demo database with 8 endmembers;
 
----- Optimization toolbox ------
qpas.*:           C-based quadratic optimization codes;
qradprogC:        A full version of the aboved optimization toolbox,  attached here for respecting its completeness. 