1

   Abaqus 6.9-EF1                                  Date 13-May-2015   Time 17:15:07
   For use at UNIVERSITE DE SAVOIE under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                1.000E-02
          AND A TOTAL TIME PERIOD OF                            1.00    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-08
          THE MAXIMUM TIME INCREMENT ALLOWED IS                1.000E-02

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM CUT-BACKS ALLOWED IN AN INCREMENT                                 5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE    6
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     AUTOMATIC TOLERANCES FOR OVERCLOSURE AND SEPARATION 
     PRESSURE ARE SUPPRESSED
     GLOBAL STABILIZATION CONTROL IS NOT USED
     FRICTION IS INCLUDED IN INCREMENT THAT THE CONTACT POINT CLOSES

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                 1000000

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.200    

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     EQUATIONS ARE BEING REORDERED TO MINIMIZE WAVEFRONT

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.492       TIME AVG. FORCE       0.492    
 LARGEST RESIDUAL FORCE              1.01       AT NODE        130   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.501       TIME AVG. FORCE       0.501    
 LARGEST RESIDUAL FORCE             1.508E-02   AT NODE        130   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.315E-05   AT NODE         12   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      0.501       TIME AVG. FORCE       0.501    
 LARGEST RESIDUAL FORCE            -1.342E-05   AT NODE         57   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.512E-06   AT NODE         21   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  1.000E-02
 STEP TIME COMPLETED       1.000E-02,  TOTAL TIME COMPLETED        1.000E-02


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.998       TIME AVG. FORCE       0.750    
 LARGEST RESIDUAL FORCE             0.277       AT NODE        121   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -3.357E-04   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.996       TIME AVG. FORCE       0.749    
 LARGEST RESIDUAL FORCE            -1.148E-02   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -9.704E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      0.997       TIME AVG. FORCE       0.749    
 LARGEST RESIDUAL FORCE            -4.067E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE         31   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -5.381E-06   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  2.000E-02
 STEP TIME COMPLETED       2.000E-02,  TOTAL TIME COMPLETED        2.000E-02


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.34       TIME AVG. FORCE       0.947    
 LARGEST RESIDUAL FORCE            -0.119       AT NODE        176   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE         30   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -3.841E-04   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.34       TIME AVG. FORCE       0.947    
 LARGEST RESIDUAL FORCE            -1.489E-03   AT NODE        176   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.207E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -6.709E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.431E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  3.000E-02
 STEP TIME COMPLETED       3.000E-02,  TOTAL TIME COMPLETED        3.000E-02


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.54       TIME AVG. FORCE        1.09    
 LARGEST RESIDUAL FORCE            -9.252E-02   AT NODE        197   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.585E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -3.783E-04   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.54       TIME AVG. FORCE        1.09    
 LARGEST RESIDUAL FORCE            -4.214E-04   AT NODE        197   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.613E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -2.777E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.265E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  4.000E-02
 STEP TIME COMPLETED       4.000E-02,  TOTAL TIME COMPLETED        4.000E-02


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.65       TIME AVG. FORCE        1.21    
 LARGEST RESIDUAL FORCE            -3.494E-02   AT NODE        208   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.667E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.329E-04   AT NODE         54   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.66       TIME AVG. FORCE        1.21    
 LARGEST RESIDUAL FORCE             4.224E-05   AT NODE        208   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.668E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.399E-06   AT NODE         31   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  5.000E-02
 STEP TIME COMPLETED       5.000E-02,  TOTAL TIME COMPLETED        5.000E-02


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.72       TIME AVG. FORCE        1.29    
 LARGEST RESIDUAL FORCE             5.535E-02   AT NODE         79   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.658E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.087E-04   AT NODE         19   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.72       TIME AVG. FORCE        1.29    
 LARGEST RESIDUAL FORCE            -5.739E-05   AT NODE         79   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.657E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.399E-06   AT NODE         79   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  6.000E-02
 STEP TIME COMPLETED       6.000E-02,  TOTAL TIME COMPLETED        6.000E-02


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.76       TIME AVG. FORCE        1.36    
 LARGEST RESIDUAL FORCE             7.851E-03   AT NODE        155   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.643E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.211E-04   AT NODE         32   DOF  1
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.76       TIME AVG. FORCE        1.36    
 LARGEST RESIDUAL FORCE            -2.069E-05   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.643E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.353E-06   AT NODE         31   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  7.000E-02
 STEP TIME COMPLETED       7.000E-02,  TOTAL TIME COMPLETED        7.000E-02


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.79       TIME AVG. FORCE        1.41    
 LARGEST RESIDUAL FORCE            -4.234E-03   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.657E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.062E-05   AT NODE        192   DOF  2
   INSTANCE: ISAMPLE                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.397E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  8.000E-02
 STEP TIME COMPLETED       8.000E-02,  TOTAL TIME COMPLETED        8.000E-02


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.81       TIME AVG. FORCE        1.46    
 LARGEST RESIDUAL FORCE            -1.436E-03   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.656E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -2.492E-05   AT NODE        145   DOF  3
   INSTANCE: ISAMPLE                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.193E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  9.000E-02
 STEP TIME COMPLETED       9.000E-02,  TOTAL TIME COMPLETED        9.000E-02


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.82       TIME AVG. FORCE        1.49    
 LARGEST RESIDUAL FORCE             8.118E-03   AT NODE        158   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.655E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.199E-05   AT NODE        126   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.82       TIME AVG. FORCE        1.49    
 LARGEST RESIDUAL FORCE             1.229E-05   AT NODE        152   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.655E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.230E-07   AT NODE        158   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED        0.100    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.83       TIME AVG. FORCE        1.52    
 LARGEST RESIDUAL FORCE            -9.563E-04   AT NODE        194   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.644E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.156E-05   AT NODE        188   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.110    
 STEP TIME COMPLETED       0.110    ,  TOTAL TIME COMPLETED        0.110    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.84       TIME AVG. FORCE        1.55    
 LARGEST RESIDUAL FORCE             2.793E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.631E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.285E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.120    
 STEP TIME COMPLETED       0.120    ,  TOTAL TIME COMPLETED        0.120    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.85       TIME AVG. FORCE        1.57    
 LARGEST RESIDUAL FORCE             2.136E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.623E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.905E-06   AT NODE        210   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.130    
 STEP TIME COMPLETED       0.130    ,  TOTAL TIME COMPLETED        0.130    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.86       TIME AVG. FORCE        1.59    
 LARGEST RESIDUAL FORCE            -1.491E-03   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.616E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.323E-06   AT NODE        209   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.140    
 STEP TIME COMPLETED       0.140    ,  TOTAL TIME COMPLETED        0.140    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.86       TIME AVG. FORCE        1.61    
 LARGEST RESIDUAL FORCE            -3.325E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.607E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -9.831E-06   AT NODE        138   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.150    
 STEP TIME COMPLETED       0.150    ,  TOTAL TIME COMPLETED        0.150    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.63    
 LARGEST RESIDUAL FORCE            -1.355E-04   AT NODE         70   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.598E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.037E-05   AT NODE         69   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.160    
 STEP TIME COMPLETED       0.160    ,  TOTAL TIME COMPLETED        0.160    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.64    
 LARGEST RESIDUAL FORCE            -1.154E-04   AT NODE        106   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.587E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.127E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.170    
 STEP TIME COMPLETED       0.170    ,  TOTAL TIME COMPLETED        0.170    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.65    
 LARGEST RESIDUAL FORCE            -9.997E-05   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.575E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.161E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.180    
 STEP TIME COMPLETED       0.180    ,  TOTAL TIME COMPLETED        0.180    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.67    
 LARGEST RESIDUAL FORCE            -9.039E-05   AT NODE         64   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.566E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.361E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.190    
 STEP TIME COMPLETED       0.190    ,  TOTAL TIME COMPLETED        0.190    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.68    
 LARGEST RESIDUAL FORCE            -1.018E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.557E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.276E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED        0.200    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.69    
 LARGEST RESIDUAL FORCE            -3.451E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.543E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.320E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.210    
 STEP TIME COMPLETED       0.210    ,  TOTAL TIME COMPLETED        0.210    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.70    
 LARGEST RESIDUAL FORCE            -1.613E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.525E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.883E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.169E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.220    
 STEP TIME COMPLETED       0.220    ,  TOTAL TIME COMPLETED        0.220    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.70    
 LARGEST RESIDUAL FORCE            -1.386E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.510E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.415E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.230    
 STEP TIME COMPLETED       0.230    ,  TOTAL TIME COMPLETED        0.230    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.71    
 LARGEST RESIDUAL FORCE            -1.390E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.499E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.216E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.240    
 STEP TIME COMPLETED       0.240    ,  TOTAL TIME COMPLETED        0.240    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.72    
 LARGEST RESIDUAL FORCE            -1.401E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.487E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.176E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.250    
 STEP TIME COMPLETED       0.250    ,  TOTAL TIME COMPLETED        0.250    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.72    
 LARGEST RESIDUAL FORCE            -1.207E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.476E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.139E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.260    
 STEP TIME COMPLETED       0.260    ,  TOTAL TIME COMPLETED        0.260    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.73    
 LARGEST RESIDUAL FORCE            -1.510E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.465E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.136E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.270    
 STEP TIME COMPLETED       0.270    ,  TOTAL TIME COMPLETED        0.270    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.74    
 LARGEST RESIDUAL FORCE            -1.409E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.455E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.093E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       0.280    ,  TOTAL TIME COMPLETED        0.280    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.74    
 LARGEST RESIDUAL FORCE            -1.471E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.444E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.066E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.290    
 STEP TIME COMPLETED       0.290    ,  TOTAL TIME COMPLETED        0.290    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.75    
 LARGEST RESIDUAL FORCE            -1.417E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.434E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.036E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED       0.300    ,  TOTAL TIME COMPLETED        0.300    


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.75    
 LARGEST RESIDUAL FORCE            -1.350E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.424E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.009E-05   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.310    
 STEP TIME COMPLETED       0.310    ,  TOTAL TIME COMPLETED        0.310    


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.76    
 LARGEST RESIDUAL FORCE            -1.415E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.415E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.614E-06   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.320    
 STEP TIME COMPLETED       0.320    ,  TOTAL TIME COMPLETED        0.320    


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.76    
 LARGEST RESIDUAL FORCE            -1.402E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.405E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.431E-06   AT NODE        212   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.330    
 STEP TIME COMPLETED       0.330    ,  TOTAL TIME COMPLETED        0.330    


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.76    
 LARGEST RESIDUAL FORCE            -1.391E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.396E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.439E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.340    
 STEP TIME COMPLETED       0.340    ,  TOTAL TIME COMPLETED        0.340    


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.77    
 LARGEST RESIDUAL FORCE            -1.411E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.387E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.162E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.350    
 STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED        0.350    


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.77    
 LARGEST RESIDUAL FORCE            -1.424E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.378E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.045E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.360    
 STEP TIME COMPLETED       0.360    ,  TOTAL TIME COMPLETED        0.360    


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.77    
 LARGEST RESIDUAL FORCE            -1.397E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.369E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.894E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.370    
 STEP TIME COMPLETED       0.370    ,  TOTAL TIME COMPLETED        0.370    


  INCREMENT    38 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.78    
 LARGEST RESIDUAL FORCE            -1.412E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.360E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.821E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.380    
 STEP TIME COMPLETED       0.380    ,  TOTAL TIME COMPLETED        0.380    


  INCREMENT    39 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.78    
 LARGEST RESIDUAL FORCE            -1.593E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.350E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.706E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.390    
 STEP TIME COMPLETED       0.390    ,  TOTAL TIME COMPLETED        0.390    


  INCREMENT    40 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.78    
 LARGEST RESIDUAL FORCE            -1.352E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.339E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.080E-05   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.400    
 STEP TIME COMPLETED       0.400    ,  TOTAL TIME COMPLETED        0.400    


  INCREMENT    41 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.78    
 LARGEST RESIDUAL FORCE            -1.344E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.330E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.453E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.410    
 STEP TIME COMPLETED       0.410    ,  TOTAL TIME COMPLETED        0.410    


  INCREMENT    42 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.79    
 LARGEST RESIDUAL FORCE            -1.428E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.321E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.336E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.420    
 STEP TIME COMPLETED       0.420    ,  TOTAL TIME COMPLETED        0.420    


  INCREMENT    43 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.79    
 LARGEST RESIDUAL FORCE            -1.358E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.312E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.736E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.430    
 STEP TIME COMPLETED       0.430    ,  TOTAL TIME COMPLETED        0.430    


  INCREMENT    44 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.79    
 LARGEST RESIDUAL FORCE            -1.355E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.303E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.563E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.440    
 STEP TIME COMPLETED       0.440    ,  TOTAL TIME COMPLETED        0.440    


  INCREMENT    45 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.79    
 LARGEST RESIDUAL FORCE            -1.438E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.295E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.432E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.450    
 STEP TIME COMPLETED       0.450    ,  TOTAL TIME COMPLETED        0.450    


  INCREMENT    46 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.80    
 LARGEST RESIDUAL FORCE            -1.359E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.287E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.290E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.460    
 STEP TIME COMPLETED       0.460    ,  TOTAL TIME COMPLETED        0.460    


  INCREMENT    47 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.80    
 LARGEST RESIDUAL FORCE            -1.198E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.278E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.238E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.470    
 STEP TIME COMPLETED       0.470    ,  TOTAL TIME COMPLETED        0.470    


  INCREMENT    48 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.80    
 LARGEST RESIDUAL FORCE            -1.320E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.270E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.092E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.480    
 STEP TIME COMPLETED       0.480    ,  TOTAL TIME COMPLETED        0.480    


  INCREMENT    49 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.80    
 LARGEST RESIDUAL FORCE            -1.335E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.262E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.977E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.490    
 STEP TIME COMPLETED       0.490    ,  TOTAL TIME COMPLETED        0.490    


  INCREMENT    50 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.80    
 LARGEST RESIDUAL FORCE            -1.327E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.254E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.871E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED       0.500    ,  TOTAL TIME COMPLETED        0.500    


  INCREMENT    51 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.317E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.246E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.302E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.510    
 STEP TIME COMPLETED       0.510    ,  TOTAL TIME COMPLETED        0.510    


  INCREMENT    52 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.305E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.237E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.005E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.520    
 STEP TIME COMPLETED       0.520    ,  TOTAL TIME COMPLETED        0.520    


  INCREMENT    53 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.305E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.229E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.949E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.530    
 STEP TIME COMPLETED       0.530    ,  TOTAL TIME COMPLETED        0.530    


  INCREMENT    54 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.298E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.221E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.732E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.540    
 STEP TIME COMPLETED       0.540    ,  TOTAL TIME COMPLETED        0.540    


  INCREMENT    55 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.292E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.214E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.563E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.550    
 STEP TIME COMPLETED       0.550    ,  TOTAL TIME COMPLETED        0.550    


  INCREMENT    56 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.286E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.        -1.206E-03   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.429E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.560    
 STEP TIME COMPLETED       0.560    ,  TOTAL TIME COMPLETED        0.560    


  INCREMENT    57 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.81    
 LARGEST RESIDUAL FORCE            -1.280E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.315E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.570    
 STEP TIME COMPLETED       0.570    ,  TOTAL TIME COMPLETED        0.570    


  INCREMENT    58 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.228E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.217E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.580    
 STEP TIME COMPLETED       0.580    ,  TOTAL TIME COMPLETED        0.580    


  INCREMENT    59 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.213E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.126E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.590    
 STEP TIME COMPLETED       0.590    ,  TOTAL TIME COMPLETED        0.590    


  INCREMENT    60 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.259E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        7.035E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED       0.600    ,  TOTAL TIME COMPLETED        0.600    


  INCREMENT    61 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.226E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.958E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.610    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    62 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.168E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.916E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.620    
 STEP TIME COMPLETED       0.620    ,  TOTAL TIME COMPLETED        0.620    


  INCREMENT    63 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.253E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.799E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.630    
 STEP TIME COMPLETED       0.630    ,  TOTAL TIME COMPLETED        0.630    


  INCREMENT    64 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.250E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.713E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.640    
 STEP TIME COMPLETED       0.640    ,  TOTAL TIME COMPLETED        0.640    


  INCREMENT    65 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.89       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.342E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.519E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.650    
 STEP TIME COMPLETED       0.650    ,  TOTAL TIME COMPLETED        0.650    


  INCREMENT    66 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.259E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.534E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.660    
 STEP TIME COMPLETED       0.660    ,  TOTAL TIME COMPLETED        0.660    


  INCREMENT    67 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.82    
 LARGEST RESIDUAL FORCE            -1.251E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.464E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.670    
 STEP TIME COMPLETED       0.670    ,  TOTAL TIME COMPLETED        0.670    


  INCREMENT    68 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.244E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.393E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.680    
 STEP TIME COMPLETED       0.680    ,  TOTAL TIME COMPLETED        0.680    


  INCREMENT    69 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.238E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.321E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.690    
 STEP TIME COMPLETED       0.690    ,  TOTAL TIME COMPLETED        0.690    


  INCREMENT    70 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.233E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.250E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.700    
 STEP TIME COMPLETED       0.700    ,  TOTAL TIME COMPLETED        0.700    


  INCREMENT    71 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.193E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.186E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.710    
 STEP TIME COMPLETED       0.710    ,  TOTAL TIME COMPLETED        0.710    


  INCREMENT    72 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.140E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.133E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.720    
 STEP TIME COMPLETED       0.720    ,  TOTAL TIME COMPLETED        0.720    


  INCREMENT    73 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.220E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.056E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.730    
 STEP TIME COMPLETED       0.730    ,  TOTAL TIME COMPLETED        0.730    


  INCREMENT    74 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.216E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.987E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.740    
 STEP TIME COMPLETED       0.740    ,  TOTAL TIME COMPLETED        0.740    


  INCREMENT    75 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.212E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.919E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.750    
 STEP TIME COMPLETED       0.750    ,  TOTAL TIME COMPLETED        0.750    


  INCREMENT    76 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.207E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.853E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.760    
 STEP TIME COMPLETED       0.760    ,  TOTAL TIME COMPLETED        0.760    


  INCREMENT    77 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.202E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.788E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.770    
 STEP TIME COMPLETED       0.770    ,  TOTAL TIME COMPLETED        0.770    


  INCREMENT    78 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.198E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.724E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.780    
 STEP TIME COMPLETED       0.780    ,  TOTAL TIME COMPLETED        0.780    


  INCREMENT    79 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.193E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.662E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.790    
 STEP TIME COMPLETED       0.790    ,  TOTAL TIME COMPLETED        0.790    


  INCREMENT    80 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.188E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.600E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.800    
 STEP TIME COMPLETED       0.800    ,  TOTAL TIME COMPLETED        0.800    


  INCREMENT    81 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.83    
 LARGEST RESIDUAL FORCE            -1.184E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.540E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.810    
 STEP TIME COMPLETED       0.810    ,  TOTAL TIME COMPLETED        0.810    


  INCREMENT    82 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.179E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.480E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.820    
 STEP TIME COMPLETED       0.820    ,  TOTAL TIME COMPLETED        0.820    


  INCREMENT    83 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.175E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.422E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.830    
 STEP TIME COMPLETED       0.830    ,  TOTAL TIME COMPLETED        0.830    


  INCREMENT    84 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.161E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.392E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.840    
 STEP TIME COMPLETED       0.840    ,  TOTAL TIME COMPLETED        0.840    


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.097E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.567E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.850    
 STEP TIME COMPLETED       0.850    ,  TOTAL TIME COMPLETED        0.850    


  INCREMENT    86 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.167E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.271E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.860    
 STEP TIME COMPLETED       0.860    ,  TOTAL TIME COMPLETED        0.860    


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.163E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.216E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.870    
 STEP TIME COMPLETED       0.870    ,  TOTAL TIME COMPLETED        0.870    


  INCREMENT    88 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.158E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.163E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.880    
 STEP TIME COMPLETED       0.880    ,  TOTAL TIME COMPLETED        0.880    


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.153E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.111E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.890    
 STEP TIME COMPLETED       0.890    ,  TOTAL TIME COMPLETED        0.890    


  INCREMENT    90 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.148E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.059E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.900    
 STEP TIME COMPLETED       0.900    ,  TOTAL TIME COMPLETED        0.900    


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.144E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.007E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.910    
 STEP TIME COMPLETED       0.910    ,  TOTAL TIME COMPLETED        0.910    


  INCREMENT    92 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.88       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.139E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.956E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.920    
 STEP TIME COMPLETED       0.920    ,  TOTAL TIME COMPLETED        0.920    


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.135E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.906E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.930    
 STEP TIME COMPLETED       0.930    ,  TOTAL TIME COMPLETED        0.930    


  INCREMENT    94 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.130E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.856E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.940    
 STEP TIME COMPLETED       0.940    ,  TOTAL TIME COMPLETED        0.940    


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.126E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.808E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.950    
 STEP TIME COMPLETED       0.950    ,  TOTAL TIME COMPLETED        0.950    


  INCREMENT    96 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.121E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.760E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.960    
 STEP TIME COMPLETED       0.960    ,  TOTAL TIME COMPLETED        0.960    


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.117E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.712E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.970    
 STEP TIME COMPLETED       0.970    ,  TOTAL TIME COMPLETED        0.970    


  INCREMENT    98 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.113E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.666E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.980    
 STEP TIME COMPLETED       0.980    ,  TOTAL TIME COMPLETED        0.980    


  INCREMENT    99 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.109E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.620E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.990    
 STEP TIME COMPLETED       0.990    ,  TOTAL TIME COMPLETED        0.990    


  INCREMENT   100 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       1.87       TIME AVG. FORCE        1.84    
 LARGEST RESIDUAL FORCE            -1.105E-04   AT NODE          6   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE        174   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.575E-06   AT NODE        211   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF        100  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                     110  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     110  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     110  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       0  WARNING MESSAGES DURING USER INPUT PROCESSING
                       0  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   5.1700    
       SYSTEM TIME (SEC)    =  0.14000    
       TOTAL CPU TIME (SEC) =   5.3100    
       WALLCLOCK TIME (SEC) =          5
