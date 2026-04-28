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
      SummarizedExperiment::assay(res)
    Output
                    graphID proteinNr   error_min    RiLog RiLog_min RiLog_max Ci
      prot_4;prot_5      NA         1 0.001632667 1.056333        NA        NA  1
                    Ci_min Ci_max case
      prot_4;prot_5     NA     NA    3

---

    Code
      SummarizedExperiment::assay(res2)
    Output
             graphID proteinNr error_min  RiLog RiLog_min RiLog_max Ci Ci_min Ci_max
      prot_1      NA         1 0.0038845 0.9635        NA        NA NA    0.1    0.9
      prot_2      NA         2 0.0038845     NA 0.9851504  1.147567 NA    0.1    0.9
             case
      prot_1    1
      prot_2    2

---

    Code
      SummarizedExperiment::assay(res3)
    Output
                    graphID proteinNr  error_min     RiLog RiLog_min RiLog_max  Ci
      prot_3             NA         1 0.00663525 1.0260321        NA        NA 0.1
      prot_4;prot_5      NA         2 0.00663525 0.9835544        NA        NA 0.9
                    Ci_min Ci_max case
      prot_3            NA     NA    3
      prot_4;prot_5     NA     NA    3

---

    Code
      SummarizedExperiment::assay(res4)
    Output
                    graphID proteinNr  error_min    RiLog RiLog_min RiLog_max Ci
      prot_3             NA         1 0.00663525 1.026032        NA        NA NA
      prot_4;prot_5      NA         2 0.00663525 0.983554        NA        NA NA
                    Ci_min Ci_max case
      prot_3           0.1    0.1    4
      prot_4;prot_5    0.9    0.9    4

---

    Code
      SummarizedExperiment::assay(res5)
    Output
                    graphID proteinNr   error_min     RiLog RiLog_min RiLog_max  Ci
      prot_3             NA         1 0.006932845 0.9779656        NA        NA 0.9
      prot_4;prot_5      NA         2 0.006635250 0.9835544        NA        NA 0.9
                    Ci_min Ci_max case
      prot_3            NA     NA    3
      prot_4;prot_5     NA     NA    3

