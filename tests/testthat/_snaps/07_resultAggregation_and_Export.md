# result aggregation

    Code
      X
    Output
      class: SummarizedExperiment 
      dim: 4 10 
      metadata(0):
      assays(3): sample1_sample2 sample1_sample3 sample2_sample3
      rownames(4): prot_1 prot_2 prot_3 prot_4;prot_5
      rowData names(1): accession
      colnames(10): graphID proteinNr ... Ci_max case
      colData names(1): colnames

---

    Code
      SummarizedExperiment::assays(X)$sample1_sample2
    Output
                    graphID proteinNr    error_min    RiLog RiLog_min RiLog_max Ci
      prot_1             NA         1 1.750500e-03 0.967500        NA        NA NA
      prot_2             NA         2 1.750500e-03       NA  1.045762  3.672106 NA
      prot_3             NA         1 1.232595e-32 0.964000        NA        NA  1
      prot_4;prot_5      NA         1 1.632667e-03 1.056333        NA        NA  1
                    Ci_min Ci_max case
      prot_1          0.01   0.99    1
      prot_2          0.01   0.99    2
      prot_3            NA     NA    3
      prot_4;prot_5     NA     NA    3

---

    Code
      SummarizedExperiment::assays(X)$sample1_sample3
    Output
                    graphID proteinNr   error_min     RiLog RiLog_min RiLog_max   Ci
      prot_1             NA         1 0.000264500 0.9910000        NA        NA   NA
      prot_2             NA         2 0.000264500        NA 0.9975651  1.528617   NA
      prot_3             NA         2 0.005563039 1.0854065        NA        NA 0.01
      prot_4;prot_5      NA         1 0.005563039 0.9528013        NA        NA 0.99
                    Ci_min Ci_max case
      prot_1          0.01   0.99    1
      prot_2          0.01   0.99    2
      prot_3            NA     NA    3
      prot_4;prot_5     NA     NA    3

---

    Code
      SummarizedExperiment::assays(X)$sample2_sample3
    Output
                    graphID proteinNr   error_min     RiLog RiLog_min RiLog_max   Ci
      prot_1             NA         1 0.003884500        NA 0.9634997 0.9635009   NA
      prot_2             NA         2 0.003884500        NA 0.9831944 2.2027825   NA
      prot_3             NA         1 0.006082681 1.0323458        NA        NA 0.01
      prot_4;prot_5      NA         2 0.006082681 0.9827164        NA        NA 0.99
                    Ci_min Ci_max case
      prot_1          0.01   0.99    2
      prot_2          0.01   0.99    2
      prot_3            NA     NA    3
      prot_4;prot_5     NA     NA    3

---

    Code
      SummarizedExperiment::rowData(X)
    Output
      DataFrame with 4 rows and 1 column
                        accession
                      <character>
      prot_1               prot_1
      prot_2               prot_2
      prot_3               prot_3
      prot_4;prot_5 prot_4;prot_5

---

    Code
      SummarizedExperiment::colData(X)
    Output
      DataFrame with 10 rows and 1 column
                   colnames
                <character>
      graphID       graphID
      proteinNr   proteinNr
      error_min   error_min
      RiLog           RiLog
      RiLog_min   RiLog_min
      RiLog_max   RiLog_max
      Ci                 Ci
      Ci_min         Ci_min
      Ci_max         Ci_max
      case             case

