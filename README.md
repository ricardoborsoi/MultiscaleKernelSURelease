#  A Blind Multiscale Spatial Regularization Framework for Kernel-Based Spectral Unmixing    #

This package contains the authors' implementation of the paper [1].

We consider a multiscale strategy using superpixels to introduce spatial information on the abundance maps for nonlinear spectral unmixing based on kernels. We split the unmixing problem in two, at a "coarse" and at a "fine" spatial scale, and formulate them using nonconvex quadratic equality constraints, which allows for the estimation of the regularization parameters blindly from the data. To solve the cost functions, we dualize the optimization problem and apply a root-finding strategy to obtain an efficient algorithm. This solution is optimal since we are able to show that strong duality holds for these optimization problems.

The one parameter which needs to be adjusted for this algorithm is related to the modeling errors.

The code is implemented in MATLAB and includes:  
-  synthetic_examples.m      - a demo script comparing the algorithms using synthetic data  
-  demo_cuprite.m            - a demo script comparing the algorithms for the Cuprite image  
-  demo_jasperRidge.m        - a demo script comparing the algorithms for the Jasper Ridge image  
-  demo_urban.m              - a demo script comparing the algorithms for the Urban subimage  
-  ./BMUAN/                  - contains the MATLAB files associated with the BMUA-N algorithm  
-  ./other_methods/          - contains the ..... methods  
-  ./utils/                  - useful functions  
-  ./DATA/                   - images used in the examples  
-  README                    - this file  



## IMPORTANT:
If you use this software please cite the following in any resulting
publication:

    [1] A Blind Multiscale Spatial Regularization Framework for Kernel-Based Spectral Unmixing 
        R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
        IEEE Transactions on Image Processing, 2020.



## INSTALLING & RUNNING:
Just start MATLAB and run one of the demo scripts (e.g. synthetic_examples.m, etc).


## NOTES:

1.  The codes for the KHype and KHype-TV algorithms were provided by Jie Chen, available at:  
    http://www.cedric-richard.fr/Matlab/chen2013nonlinear.zip  
    http://www.cedric-richard.fr/Matlab/chen2014nonlinear2.zip

2.  The codes for the CDA-NL algorithm were provided by Abderrahim Halimi, available at:  
    https://sites.google.com/site/abderrahimhalimi/publications

3.  The codes for the NDU algorithm was provided by Rita Ammanouil, available at:  
    https://github.com/rammanouil/NonlinearUnmixingVector



