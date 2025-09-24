# test .errorEquation

    Code
      

# test .minimizeSquaredError

    Code
      res
    Output
      $Ri
      [1] 0.4787964 1.2000000
      
      $Ci
      [1] 0.8319426 0.1680574
      
      $RES
      $RES$res_Mat
                [,1]      [,2]
      [1,] 0.3983311 0.2016689
      [2,] 0.0000000 1.2000000
      
      $RES$res_equ
      [1] 2.464725e-08 2.212424e-08
      
      $RES$res_squ_err
      [1] 1.096969e-15
      
      $RES$W
                [,1]      [,2]
      [1,] 0.8319426 0.1680574
      [2,] 0.0000000 1.0000000
      
      
      $Tracking
        iter      squ_err         R1        R2        C1        C2
      1    0 1.644020e-01 -0.7369656 0.2630344 0.5000000 0.5000000
      2    1 1.096969e-15 -1.0625158 0.2630344 0.8319426 0.1680574
      
      $outer.iter
      [1] 2
      
      $convergence
      [1] 0
      
      attr(,"class")
      [1] "res_min_squ_error"

---

    Code
      res2
    Output
      $Ri
      [1] 1.231505e-08 1.014185e+00
      
      $Ci
      [1] 0.3 0.7
      
      $RES
      $RES$res_Mat
                   [,1]      [,2]
      [1,] 3.694514e-09 0.7099297
      [2,] 0.000000e+00 1.0141852
      
      $RES$res_equ
      [1] -0.1682362  0.1682360
      
      $RES$res_squ_err
      [1] 0.05660678
      
      $RES$W
           [,1] [,2]
      [1,]  0.3  0.7
      [2,]  0.0  1.0
      
      
      $Tracking
        iter    squ_err          R1         R2  C1  C2
      1    0 0.28156634  -0.7369656 0.26303441 0.3 0.7
      2    1 0.05660678 -26.2750027 0.02032116 0.3 0.7
      
      $outer.iter
      [1] 2
      
      $convergence
      [1] 0
      
      attr(,"class")
      [1] "res_min_squ_error"

