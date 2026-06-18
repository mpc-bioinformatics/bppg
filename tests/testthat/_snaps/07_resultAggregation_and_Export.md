# result aggregation

    Code
      X
    Output
      class: SummarizedExperiment 
      dim: 5 10 
      metadata(0):
      assays(36): 1_2 1_3 ... 7_9 8_9
      rownames(5): sp|P07262|DHE4_YEAST sp|P09938|RIR2_YEAST
        sp|P39708|DHE5_YEAST sp|P40212|RL13B_YEAST sp|Q12690|RL13A_YEAST
      rowData names(1): accession
      colnames(10): graphID proteinNr ... Ci_max case
      colData names(1): colnames

---

    Code
      SummarizedExperiment::assays(X)$"1_2"
    Output
                            graphID proteinNr error_min       RiLog  RiLog_min
      sp|P07262|DHE4_YEAST       NA         2  1.264748 -0.04778300         NA
      sp|P09938|RIR2_YEAST       NA         1  2.011747 -0.01770424         NA
      sp|P39708|DHE5_YEAST       NA         1  1.264748  0.21612300         NA
      sp|P40212|RL13B_YEAST      NA         1  1.302607          NA  0.5632699
      sp|Q12690|RL13A_YEAST      NA         2  1.302607          NA -5.2291086
                              RiLog_max Ci Ci_min Ci_max case
      sp|P07262|DHE4_YEAST           NA NA   0.35   0.35    4
      sp|P09938|RIR2_YEAST           NA  1     NA     NA    3
      sp|P39708|DHE5_YEAST           NA NA   0.65   0.65    4
      sp|P40212|RL13B_YEAST  0.56327184 NA   0.01   0.63    5
      sp|Q12690|RL13A_YEAST -0.09643997 NA   0.37   0.99    5

---

    Code
      SummarizedExperiment::rowData(X)
    Output
      DataFrame with 5 rows and 1 column
                                        accession
                                      <character>
      sp|P07262|DHE4_YEAST   sp|P07262|DHE4_YEAST
      sp|P09938|RIR2_YEAST   sp|P09938|RIR2_YEAST
      sp|P39708|DHE5_YEAST   sp|P39708|DHE5_YEAST
      sp|P40212|RL13B_YEAST sp|P40212|RL13B_YEAST
      sp|Q12690|RL13A_YEAST sp|Q12690|RL13A_YEAST

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

