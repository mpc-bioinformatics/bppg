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
      res
    Output
            accession comparison graphID proteinNr   error_min    RiLog RiLog_min
      1 prot_4;prot_5         NA      NA         1 0.001632667 1.056333        NA
        RiLog_max Ci Ci_min Ci_max case
      1        NA  1     NA     NA    3

---

    Code
      res2
    Output
        accession comparison graphID proteinNr error_min  RiLog RiLog_min RiLog_max
      1    prot_1         NA      NA         1 0.0038845 0.9635        NA        NA
      2    prot_2         NA      NA         2 0.0038845     NA 0.9851504  1.147567
        Ci Ci_min Ci_max case
      1 NA    0.1    0.9    1
      2 NA    0.1    0.9    2

---

    Code
      res3
    Output
            accession comparison graphID proteinNr  error_min     RiLog RiLog_min
      1        prot_3         NA      NA         1 0.00663525 1.0260321        NA
      2 prot_4;prot_5         NA      NA         2 0.00663525 0.9835544        NA
        RiLog_max  Ci Ci_min Ci_max case
      1        NA 0.1     NA     NA    3
      2        NA 0.9     NA     NA    3

---

    Code
      res4
    Output
            accession comparison graphID proteinNr  error_min    RiLog RiLog_min
      1        prot_3         NA      NA         1 0.00663525 1.026032        NA
      2 prot_4;prot_5         NA      NA         2 0.00663525 0.983554        NA
        RiLog_max Ci Ci_min Ci_max case
      1        NA NA    0.1    0.1    4
      2        NA NA    0.9    0.9    4

