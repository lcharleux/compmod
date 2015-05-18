1

   Abaqus 6.12-1                                  Date 15-mai-2015   Time 14:45:08
   For use at  under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
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
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
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
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                 1000000

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     5.000E-02

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     EQUATIONS ARE BEING REORDERED TO MINIMIZE WAVEFRONT

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  6.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2700    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.073E-02   TIME AVG. FORCE       3.073E-02
 LARGEST RESIDUAL FORCE             -1.41       AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2500    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2700    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.147E-02   TIME AVG. FORCE       3.147E-02
 LARGEST RESIDUAL FORCE             1.509E-02   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -5.190E-05   AT NODE        400   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.147E-02   TIME AVG. FORCE       3.147E-02
 LARGEST RESIDUAL FORCE             3.432E-04   AT NODE        421   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.655E-05   AT NODE        421   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.4200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.4300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.147E-02   TIME AVG. FORCE       3.147E-02
 LARGEST RESIDUAL FORCE            -6.277E-06   AT NODE       1237   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.950E-06   AT NODE        421   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  1.000E-02
 STEP TIME COMPLETED       1.000E-02,  TOTAL TIME COMPLETED        1.000E-02


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.216E-02   TIME AVG. FORCE       4.682E-02
 LARGEST RESIDUAL FORCE            -9.686E-02   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.271E-04   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.215E-02   TIME AVG. FORCE       4.681E-02
 LARGEST RESIDUAL FORCE             1.705E-03   AT NODE        909   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.745E-05   AT NODE          5   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.215E-02   TIME AVG. FORCE       4.681E-02
 LARGEST RESIDUAL FORCE             4.977E-06   AT NODE        447   DOF  3
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.820E-06   AT NODE          7   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  2.000E-02
 STEP TIME COMPLETED       2.000E-02,  TOTAL TIME COMPLETED        2.000E-02


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.726E-02   TIME AVG. FORCE       6.029E-02
 LARGEST RESIDUAL FORCE            -0.149       AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.201E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.654E-04   AT NODE       8389   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.720E-02   TIME AVG. FORCE       6.027E-02
 LARGEST RESIDUAL FORCE            -1.498E-03   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.575E-05   AT NODE        211   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.720E-02   TIME AVG. FORCE       6.027E-02
 LARGEST RESIDUAL FORCE            -1.275E-06   AT NODE       8654   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -3.126E-07   AT NODE       6783   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  3.000E-02
 STEP TIME COMPLETED       3.000E-02,  TOTAL TIME COMPLETED        3.000E-02


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.104       TIME AVG. FORCE       7.128E-02
 LARGEST RESIDUAL FORCE            -7.637E-02   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.218E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.034E-04   AT NODE       8840   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.104       TIME AVG. FORCE       7.126E-02
 LARGEST RESIDUAL FORCE             1.040E-03   AT NODE       8926   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        8.943E-06   AT NODE       9175   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      0.104       TIME AVG. FORCE       7.126E-02
 LARGEST RESIDUAL FORCE             1.582E-06   AT NODE       8905   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -8.056E-07   AT NODE       8905   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  4.000E-02
 STEP TIME COMPLETED       4.000E-02,  TOTAL TIME COMPLETED        4.000E-02


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.115       TIME AVG. FORCE       7.997E-02
 LARGEST RESIDUAL FORCE            -0.109       AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.218E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -9.335E-05   AT NODE       8838   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2300    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2600    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.115       TIME AVG. FORCE       7.997E-02
 LARGEST RESIDUAL FORCE             1.750E-04   AT NODE       7679   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.218E-03   AT NODE       9240   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -4.873E-06   AT NODE       7098   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  5.000E-02
 STEP TIME COMPLETED       5.000E-02,  TOTAL TIME COMPLETED        5.000E-02


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       8.685E-02
 LARGEST RESIDUAL FORCE            -4.453E-02   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.217E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -5.590E-05   AT NODE       8925   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       8.685E-02
 LARGEST RESIDUAL FORCE             9.913E-05   AT NODE       8524   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.217E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -2.013E-06   AT NODE       4641   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  6.000E-02
 STEP TIME COMPLETED       6.000E-02,  TOTAL TIME COMPLETED        6.000E-02


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.125       TIME AVG. FORCE       9.233E-02
 LARGEST RESIDUAL FORCE            -1.057E-02   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -2.167E-05   AT NODE       9187   DOF  3
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.125       TIME AVG. FORCE       9.233E-02
 LARGEST RESIDUAL FORCE             1.238E-05   AT NODE       5766   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.299E-06   AT NODE       7466   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  7.000E-02
 STEP TIME COMPLETED       7.000E-02,  TOTAL TIME COMPLETED        7.000E-02


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.128       TIME AVG. FORCE       9.675E-02
 LARGEST RESIDUAL FORCE            -5.116E-03   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.073E-05   AT NODE        195   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.128       TIME AVG. FORCE       9.675E-02
 LARGEST RESIDUAL FORCE            -2.982E-05   AT NODE        561   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.369E-06   AT NODE       7032   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  8.000E-02
 STEP TIME COMPLETED       8.000E-02,  TOTAL TIME COMPLETED        8.000E-02


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2300    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2500    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.129       TIME AVG. FORCE       0.100    
 LARGEST RESIDUAL FORCE            -4.876E-03   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.385E-05   AT NODE       1927   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.129       TIME AVG. FORCE       0.100    
 LARGEST RESIDUAL FORCE            -5.035E-06   AT NODE       5123   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -8.376E-07   AT NODE       5544   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  9.000E-02
 STEP TIME COMPLETED       9.000E-02,  TOTAL TIME COMPLETED        9.000E-02


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.130       TIME AVG. FORCE       0.103    
 LARGEST RESIDUAL FORCE            -2.373E-03   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -1.136E-05   AT NODE        115   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.130       TIME AVG. FORCE       0.103    
 LARGEST RESIDUAL FORCE            -3.635E-06   AT NODE       5695   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -7.064E-07   AT NODE         45   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED        0.100    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.131       TIME AVG. FORCE       0.106    
 LARGEST RESIDUAL FORCE             7.432E-04   AT NODE       2462   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        9.016E-06   AT NODE       9038   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.131       TIME AVG. FORCE       0.106    
 LARGEST RESIDUAL FORCE            -2.809E-06   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.219E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.784E-07   AT NODE       2462   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.110    
 STEP TIME COMPLETED       0.110    ,  TOTAL TIME COMPLETED        0.110    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.108    
 LARGEST RESIDUAL FORCE            -1.652E-03   AT NODE       1585   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.217E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.511E-06   AT NODE        389   DOF  2
   INSTANCE: ISAMPLE                                                                         
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.108    
 LARGEST RESIDUAL FORCE             3.394E-06   AT NODE       1585   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.217E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -4.436E-07   AT NODE       1585   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.120    
 STEP TIME COMPLETED       0.120    ,  TOTAL TIME COMPLETED        0.120    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.110    
 LARGEST RESIDUAL FORCE             4.333E-04   AT NODE       1605   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.215E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        5.534E-06   AT NODE       6619   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.130    
 STEP TIME COMPLETED       0.130    ,  TOTAL TIME COMPLETED        0.130    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.112    
 LARGEST RESIDUAL FORCE             1.121E-04   AT NODE       7655   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.214E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -5.461E-06   AT NODE       6148   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.140    
 STEP TIME COMPLETED       0.140    ,  TOTAL TIME COMPLETED        0.140    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2400    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2500    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.113    
 LARGEST RESIDUAL FORCE             1.179E-04   AT NODE       8888   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.213E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        4.368E-06   AT NODE       8995   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.150    
 STEP TIME COMPLETED       0.150    ,  TOTAL TIME COMPLETED        0.150    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.114    
 LARGEST RESIDUAL FORCE            -1.573E-04   AT NODE       7828   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.213E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        6.857E-06   AT NODE       9009   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.160    
 STEP TIME COMPLETED       0.160    ,  TOTAL TIME COMPLETED        0.160    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.115    
 LARGEST RESIDUAL FORCE             3.940E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.212E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        3.059E-06   AT NODE       9255   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.170    
 STEP TIME COMPLETED       0.170    ,  TOTAL TIME COMPLETED        0.170    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.116    
 LARGEST RESIDUAL FORCE             3.108E-04   AT NODE       8479   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.212E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        3.528E-06   AT NODE        302   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.180    
 STEP TIME COMPLETED       0.180    ,  TOTAL TIME COMPLETED        0.180    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.117    
 LARGEST RESIDUAL FORCE            -5.300E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.211E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.919E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.190    
 STEP TIME COMPLETED       0.190    ,  TOTAL TIME COMPLETED        0.190    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.118    
 LARGEST RESIDUAL FORCE            -6.564E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.211E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.816E-06   AT NODE       9065   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED        0.200    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.119    
 LARGEST RESIDUAL FORCE            -5.166E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.211E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.590E-06   AT NODE         36   DOF  2
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.210    
 STEP TIME COMPLETED       0.210    ,  TOTAL TIME COMPLETED        0.210    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.119    
 LARGEST RESIDUAL FORCE            -7.773E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.210E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.       -2.341E-06   AT NODE          5   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.220    
 STEP TIME COMPLETED       0.220    ,  TOTAL TIME COMPLETED        0.220    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.120    
 LARGEST RESIDUAL FORCE            -8.690E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.210E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.448E-06   AT NODE       8841   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.230    
 STEP TIME COMPLETED       0.230    ,  TOTAL TIME COMPLETED        0.230    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.121    
 LARGEST RESIDUAL FORCE            -9.044E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.209E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.306E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.240    
 STEP TIME COMPLETED       0.240    ,  TOTAL TIME COMPLETED        0.240    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.121    
 LARGEST RESIDUAL FORCE            -9.872E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.209E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.280E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.250    
 STEP TIME COMPLETED       0.250    ,  TOTAL TIME COMPLETED        0.250    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.122    
 LARGEST RESIDUAL FORCE            -8.961E-05   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.209E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.280E-06   AT NODE       7119   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.260    
 STEP TIME COMPLETED       0.260    ,  TOTAL TIME COMPLETED        0.260    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.122    
 LARGEST RESIDUAL FORCE            -1.102E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.208E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.251E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.270    
 STEP TIME COMPLETED       0.270    ,  TOTAL TIME COMPLETED        0.270    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.122    
 LARGEST RESIDUAL FORCE            -1.110E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.208E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.238E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       0.280    ,  TOTAL TIME COMPLETED        0.280    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.123    
 LARGEST RESIDUAL FORCE            -1.139E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.208E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.194E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.290    
 STEP TIME COMPLETED       0.290    ,  TOTAL TIME COMPLETED        0.290    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.3800    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.4000    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.123    
 LARGEST RESIDUAL FORCE            -1.158E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.207E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.182E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED       0.300    ,  TOTAL TIME COMPLETED        0.300    


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2700    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2800    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.124    
 LARGEST RESIDUAL FORCE            -1.162E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.207E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.137E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.310    
 STEP TIME COMPLETED       0.310    ,  TOTAL TIME COMPLETED        0.310    


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.3600    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.3800    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.124    
 LARGEST RESIDUAL FORCE            -1.168E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.207E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.135E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.320    
 STEP TIME COMPLETED       0.320    ,  TOTAL TIME COMPLETED        0.320    


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.3000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.3100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.124    
 LARGEST RESIDUAL FORCE            -1.189E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.207E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.131E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.330    
 STEP TIME COMPLETED       0.330    ,  TOTAL TIME COMPLETED        0.330    


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.124    
 LARGEST RESIDUAL FORCE            -1.195E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.119E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.340    
 STEP TIME COMPLETED       0.340    ,  TOTAL TIME COMPLETED        0.340    


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.125    
 LARGEST RESIDUAL FORCE            -1.199E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.108E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.350    
 STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED        0.350    


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.125    
 LARGEST RESIDUAL FORCE            -1.154E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.206E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.096E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.360    
 STEP TIME COMPLETED       0.360    ,  TOTAL TIME COMPLETED        0.360    


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.125    
 LARGEST RESIDUAL FORCE            -1.137E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.080E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.370    
 STEP TIME COMPLETED       0.370    ,  TOTAL TIME COMPLETED        0.370    


  INCREMENT    38 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.125    
 LARGEST RESIDUAL FORCE            -1.204E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.056E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.380    
 STEP TIME COMPLETED       0.380    ,  TOTAL TIME COMPLETED        0.380    


  INCREMENT    39 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.126    
 LARGEST RESIDUAL FORCE            -1.204E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.044E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.390    
 STEP TIME COMPLETED       0.390    ,  TOTAL TIME COMPLETED        0.390    


  INCREMENT    40 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.126    
 LARGEST RESIDUAL FORCE            -1.215E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.205E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.031E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.400    
 STEP TIME COMPLETED       0.400    ,  TOTAL TIME COMPLETED        0.400    


  INCREMENT    41 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.126    
 LARGEST RESIDUAL FORCE            -1.214E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.017E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.410    
 STEP TIME COMPLETED       0.410    ,  TOTAL TIME COMPLETED        0.410    


  INCREMENT    42 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.126    
 LARGEST RESIDUAL FORCE            -1.213E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        2.004E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.420    
 STEP TIME COMPLETED       0.420    ,  TOTAL TIME COMPLETED        0.420    


  INCREMENT    43 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.126    
 LARGEST RESIDUAL FORCE            -1.211E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.991E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.430    
 STEP TIME COMPLETED       0.430    ,  TOTAL TIME COMPLETED        0.430    


  INCREMENT    44 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.210E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.204E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.971E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.440    
 STEP TIME COMPLETED       0.440    ,  TOTAL TIME COMPLETED        0.440    


  INCREMENT    45 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.208E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.964E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.450    
 STEP TIME COMPLETED       0.450    ,  TOTAL TIME COMPLETED        0.450    


  INCREMENT    46 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.205E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.960E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.460    
 STEP TIME COMPLETED       0.460    ,  TOTAL TIME COMPLETED        0.460    


  INCREMENT    47 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.203E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.948E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.470    
 STEP TIME COMPLETED       0.470    ,  TOTAL TIME COMPLETED        0.470    


  INCREMENT    48 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.200E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.203E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.934E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.480    
 STEP TIME COMPLETED       0.480    ,  TOTAL TIME COMPLETED        0.480    


  INCREMENT    49 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.197E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.936E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.490    
 STEP TIME COMPLETED       0.490    ,  TOTAL TIME COMPLETED        0.490    


  INCREMENT    50 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.194E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.941E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED       0.500    ,  TOTAL TIME COMPLETED        0.500    


  INCREMENT    51 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.127    
 LARGEST RESIDUAL FORCE            -1.191E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.904E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.510    
 STEP TIME COMPLETED       0.510    ,  TOTAL TIME COMPLETED        0.510    


  INCREMENT    52 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.187E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.898E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.520    
 STEP TIME COMPLETED       0.520    ,  TOTAL TIME COMPLETED        0.520    


  INCREMENT    53 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.183E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.201E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.888E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.530    
 STEP TIME COMPLETED       0.530    ,  TOTAL TIME COMPLETED        0.530    


  INCREMENT    54 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.179E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.201E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.874E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.540    
 STEP TIME COMPLETED       0.540    ,  TOTAL TIME COMPLETED        0.540    


  INCREMENT    55 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.176E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.201E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.861E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.550    
 STEP TIME COMPLETED       0.550    ,  TOTAL TIME COMPLETED        0.550    


  INCREMENT    56 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.173E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.201E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.852E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.560    
 STEP TIME COMPLETED       0.560    ,  TOTAL TIME COMPLETED        0.560    


  INCREMENT    57 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2500    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2500    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.169E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.842E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.570    
 STEP TIME COMPLETED       0.570    ,  TOTAL TIME COMPLETED        0.570    


  INCREMENT    58 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.165E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        402   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.831E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.580    
 STEP TIME COMPLETED       0.580    ,  TOTAL TIME COMPLETED        0.580    


  INCREMENT    59 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.162E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.821E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.590    
 STEP TIME COMPLETED       0.590    ,  TOTAL TIME COMPLETED        0.590    


  INCREMENT    60 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.158E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.811E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED       0.600    ,  TOTAL TIME COMPLETED        0.600    


  INCREMENT    61 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.154E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.801E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.610    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    62 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.128    
 LARGEST RESIDUAL FORCE            -1.150E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.791E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.620    
 STEP TIME COMPLETED       0.620    ,  TOTAL TIME COMPLETED        0.620    


  INCREMENT    63 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.146E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.781E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.630    
 STEP TIME COMPLETED       0.630    ,  TOTAL TIME COMPLETED        0.630    


  INCREMENT    64 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.142E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.772E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.640    
 STEP TIME COMPLETED       0.640    ,  TOTAL TIME COMPLETED        0.640    


  INCREMENT    65 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.138E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.774E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.650    
 STEP TIME COMPLETED       0.650    ,  TOTAL TIME COMPLETED        0.650    


  INCREMENT    66 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.134E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.772E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.660    
 STEP TIME COMPLETED       0.660    ,  TOTAL TIME COMPLETED        0.660    


  INCREMENT    67 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.130E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.757E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.670    
 STEP TIME COMPLETED       0.670    ,  TOTAL TIME COMPLETED        0.670    


  INCREMENT    68 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.126E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.746E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.680    
 STEP TIME COMPLETED       0.680    ,  TOTAL TIME COMPLETED        0.680    


  INCREMENT    69 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.4500    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.4700    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.122E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.736E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.690    
 STEP TIME COMPLETED       0.690    ,  TOTAL TIME COMPLETED        0.690    


  INCREMENT    70 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.118E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.727E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.700    
 STEP TIME COMPLETED       0.700    ,  TOTAL TIME COMPLETED        0.700    


  INCREMENT    71 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.114E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.715E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.710    
 STEP TIME COMPLETED       0.710    ,  TOTAL TIME COMPLETED        0.710    


  INCREMENT    72 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.110E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.707E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.720    
 STEP TIME COMPLETED       0.720    ,  TOTAL TIME COMPLETED        0.720    


  INCREMENT    73 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.107E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.698E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.730    
 STEP TIME COMPLETED       0.730    ,  TOTAL TIME COMPLETED        0.730    


  INCREMENT    74 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.103E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.701E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.740    
 STEP TIME COMPLETED       0.740    ,  TOTAL TIME COMPLETED        0.740    


  INCREMENT    75 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.099E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.689E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.750    
 STEP TIME COMPLETED       0.750    ,  TOTAL TIME COMPLETED        0.750    


  INCREMENT    76 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.095E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.677E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.760    
 STEP TIME COMPLETED       0.760    ,  TOTAL TIME COMPLETED        0.760    


  INCREMENT    77 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.091E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.668E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.770    
 STEP TIME COMPLETED       0.770    ,  TOTAL TIME COMPLETED        0.770    


  INCREMENT    78 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.087E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.658E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.780    
 STEP TIME COMPLETED       0.780    ,  TOTAL TIME COMPLETED        0.780    


  INCREMENT    79 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.083E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.649E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.790    
 STEP TIME COMPLETED       0.790    ,  TOTAL TIME COMPLETED        0.790    


  INCREMENT    80 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.079E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.640E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.800    
 STEP TIME COMPLETED       0.800    ,  TOTAL TIME COMPLETED        0.800    


  INCREMENT    81 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.075E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.631E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.810    
 STEP TIME COMPLETED       0.810    ,  TOTAL TIME COMPLETED        0.810    


  INCREMENT    82 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.071E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.623E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.820    
 STEP TIME COMPLETED       0.820    ,  TOTAL TIME COMPLETED        0.820    


  INCREMENT    83 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.129    
 LARGEST RESIDUAL FORCE            -1.068E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.614E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.830    
 STEP TIME COMPLETED       0.830    ,  TOTAL TIME COMPLETED        0.830    


  INCREMENT    84 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.064E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.606E-06   AT NODE       8851   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.840    
 STEP TIME COMPLETED       0.840    ,  TOTAL TIME COMPLETED        0.840    


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.060E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.598E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.850    
 STEP TIME COMPLETED       0.850    ,  TOTAL TIME COMPLETED        0.850    


  INCREMENT    86 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.056E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.592E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.860    
 STEP TIME COMPLETED       0.860    ,  TOTAL TIME COMPLETED        0.860    


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.052E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.586E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.870    
 STEP TIME COMPLETED       0.870    ,  TOTAL TIME COMPLETED        0.870    


  INCREMENT    88 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.049E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.581E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.880    
 STEP TIME COMPLETED       0.880    ,  TOTAL TIME COMPLETED        0.880    


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          3

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.045E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.575E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.890    
 STEP TIME COMPLETED       0.890    ,  TOTAL TIME COMPLETED        0.890    


  INCREMENT    90 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.041E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.570E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.900    
 STEP TIME COMPLETED       0.900    ,  TOTAL TIME COMPLETED        0.900    


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.038E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.564E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.910    
 STEP TIME COMPLETED       0.910    ,  TOTAL TIME COMPLETED        0.910    


  INCREMENT    92 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2100    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.034E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.558E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.920    
 STEP TIME COMPLETED       0.920    ,  TOTAL TIME COMPLETED        0.920    


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2300    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2400    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.030E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.553E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.930    
 STEP TIME COMPLETED       0.930    ,  TOTAL TIME COMPLETED        0.930    


  INCREMENT    94 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2300    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.027E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.547E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.940    
 STEP TIME COMPLETED       0.940    ,  TOTAL TIME COMPLETED        0.940    


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.023E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.542E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.950    
 STEP TIME COMPLETED       0.950    ,  TOTAL TIME COMPLETED        0.950    


  INCREMENT    96 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.019E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.536E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.960    
 STEP TIME COMPLETED       0.960    ,  TOTAL TIME COMPLETED        0.960    


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.1900    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.016E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.531E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.970    
 STEP TIME COMPLETED       0.970    ,  TOTAL TIME COMPLETED        0.970    


  INCREMENT    98 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2200    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.012E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.525E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.980    
 STEP TIME COMPLETED       0.980    ,  TOTAL TIME COMPLETED        0.980    


  INCREMENT    99 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.009E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.520E-06   AT NODE       8836   DOF  3
   INSTANCE: ISAMPLE                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E-02,  FRACTION OF STEP COMPLETED  0.990    
 STEP TIME COMPLETED       0.990    ,  TOTAL TIME COMPLETED        0.990    


  INCREMENT   100 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E-02

  NUMBER OF EQUATIONS =    26903     NUMBER OF FLOATING PT. OPERATIONS =  1.34E+10

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.2100    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =   2.2200    
       ELAPSED WALLCLOCK TIME (SEC) =          2

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.132       TIME AVG. FORCE       0.130    
 LARGEST RESIDUAL FORCE            -1.005E-04   AT NODE         21   DOF  1
   INSTANCE: ISAMPLE                                                                         
 LARGEST INCREMENT OF DISP.         1.200E-03   AT NODE        421   DOF  2
   INSTANCE: ISAMPLE                                                                         
 LARGEST CORRECTION TO DISP.        1.514E-06   AT NODE       8836   DOF  3
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
                     117  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     117  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     117  INVOLVE MATRIX DECOMPOSITION, INCLUDING
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
       USER TIME (SEC)      =   684.61    
       SYSTEM TIME (SEC)    =   8.5100    
       TOTAL CPU TIME (SEC) =   693.12    
       WALLCLOCK TIME (SEC) =        693
