# test .errorEquation

    Code
      e
    Output
      $res_Mat
           [,1] [,2]
      [1,] 0.50 0.00
      [2,] 0.15 0.91
      
      $res_equ
      [1] 0.2630344 0.1789701
      
      $res_squ_err
      [1] 0.1012174
      
      $W
           [,1] [,2]
      [1,]  1.0  0.0
      [2,]  0.3  0.7
      

# test automated analysis

    Code
      res3
    Output
            RiLog RiLog_min  RiLog_max Ci Ci_min Ci_max case
      1        NA -23.51333 -0.9068906 NA    0.5    0.9    5
      2 0.2630344        NA         NA NA    0.1    0.5    4

# test automated analysis with using results from other proteins

    Code
      res4
    Output
            RiLog RiLog_min  RiLog_max Ci Ci_min Ci_max case
      1        NA -23.51906 -0.9068906 NA    0.5    0.9    5
      2 0.2630344        NA         NA NA    0.1    0.5    4

