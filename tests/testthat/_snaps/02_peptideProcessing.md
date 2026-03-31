# test aggregateReplicates

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(0):
      assays(1): intensities
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): 1 2 3
      colData names(1): group

---

    Code
      SummarizedExperiment::assays(D1)$intensities
    Output
                    1        2        3
      pep_1        NA 17.44052       NA
      pep_2        NA 22.33180 18.38504
      pep_3  21.03811 16.82917       NA
      pep_4  18.20572       NA       NA
      pep_5  20.97100 15.85685 21.47879
      pep_6        NA 22.93523       NA
      pep_7  17.86183 20.95111       NA
      pep_8  18.51944 19.42229 17.56504
      pep_9  18.78308       NA 21.13843
      pep_10 20.31003 19.06793       NA

---

    Code
      as.data.frame(tail(SummarizedExperiment::rowData(D1), n = 1000))
    Output
             Sequence
      pep_1     pep_1
      pep_2     pep_2
      pep_3     pep_3
      pep_4     pep_4
      pep_5     pep_5
      pep_6     pep_6
      pep_7     pep_7
      pep_8     pep_8
      pep_9     pep_9
      pep_10   pep_10

---

    Code
      as.data.frame(tail(SummarizedExperiment::colData(D1), n = 1000))
    Output
        group
      1     1
      2     2
      3     3

---

    Code
      D2
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(0):
      assays(1): intensities
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): 1 2 3
      colData names(1): group

---

    Code
      SummarizedExperiment::assays(D2)$intensities
    Output
                    1        2        3
      pep_1  17.60185 17.18089 20.85566
      pep_2  19.96728 23.27206 18.90382
      pep_3  21.87331 16.56922 20.97479
      pep_4  18.29231 22.74236 19.09793
      pep_5  22.63031 15.89846 20.83005
      pep_6        NA 23.13152 23.35006
      pep_7  15.72499 21.84061 16.69903
      pep_8  19.18439 18.87477 17.18942
      pep_9  18.79415       NA 21.00017
      pep_10 19.92332 20.86509 19.89806

---

    Code
      as.data.frame(tail(SummarizedExperiment::rowData(D2), n = 1000))
    Output
             Sequence
      pep_1     pep_1
      pep_2     pep_2
      pep_3     pep_3
      pep_4     pep_4
      pep_5     pep_5
      pep_6     pep_6
      pep_7     pep_7
      pep_8     pep_8
      pep_9     pep_9
      pep_10   pep_10

---

    Code
      D1
    Output

      class: SummarizedExperiment 
      dim: 10 3 
      metadata(1): imputed
      assays(1): intensities
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): 1 2 3
      colData names(1): group
---

    Code
      D2
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(1): imputed
      assays(1): intensities
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): 1 2 3
      colData names(1): group

---

    Code
      as.data.frame(tail(SummarizedExperiment::colData(D2), n = 1000))
    Output
        group
      1     1
      2     2
      3     3

# test calculatePeptideRatios

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(0):
      assays(1): logRatios
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): logRatio_sample1_sample2 logRatio_sample1_sample3
        logRatio_sample2_sample3
      colData names(1): comparison

---

    Code
      SummarizedExperiment::assays(D3)$intensities
    Output
                     1        2         3
      pep_1   8.251585 17.44052  9.232288
      pep_2   7.918788 22.33180 18.385040
      pep_3  21.038106 16.82917  9.693175
      pep_4  18.205715 10.29570  9.282804
      pep_5  20.970997 15.85685 21.478791
      pep_6   9.508241 22.93523 11.012259
      pep_7  17.861831 20.95111  7.529904
      pep_8  18.519439 19.42229 17.565035
      pep_9  18.783078 10.48945 21.138431
      pep_10 20.310029 19.06793  9.325311

---

    Code
      D3
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(1): imputed
      assays(2): intensities maskImputation
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): 1 2 3
      colData names(1): group

# test calculatePeptideRatios

    Code
      SummarizedExperiment::assays(D1)$logRatios
    Output
             logRatio_sample1_sample2 logRatio_sample1_sample3
      pep_1                        NA                       NA
      pep_2                        NA               0.30166498
      pep_3                        NA                       NA
      pep_4               -0.52979568              -0.39790351
      pep_5                0.52218804               0.30487036
      pep_6                0.02686824              -0.18866120
      pep_7               -0.58551545                       NA
      pep_8                0.11089975                       NA
      pep_9                        NA              -0.03532842
      pep_10               0.39241543                       NA
             logRatio_sample2_sample3
      pep_1                        NA
      pep_2                        NA
      pep_3                        NA
      pep_4                 0.1318922
      pep_5                -0.2173177
      pep_6                -0.2155294
      pep_7                        NA
      pep_8                        NA
      pep_9                        NA
      pep_10                       NA

---

    Code
      SummarizedExperiment::rowData(D1)
    Output
      DataFrame with 10 rows and 1 column
                Sequence
             <character>
      pep_1        pep_1
      pep_2        pep_2
      pep_3        pep_3
      pep_4        pep_4
      pep_5        pep_5
      pep_6        pep_6
      pep_7        pep_7
      pep_8        pep_8
      pep_9        pep_9
      pep_10      pep_10

---

    Code
      SummarizedExperiment::colData(D1)
    Output
      DataFrame with 3 rows and 1 column
                                           comparison
                                          <character>
      logRatio_sample1_sample2 logRatio_sample1_sam..
      logRatio_sample1_sample3 logRatio_sample1_sam..
      logRatio_sample2_sample3 logRatio_sample2_sam..

# normalize peptide data

    Code
      D_norm_loess
    Output
      class: SummarizedExperiment 
      dim: 7944 27 
      metadata(0):
      assays(1): intensities_norm
      rownames(7944): AAAAQDEITGDGTTTVVCLVGELLR AAADALSDLEIK ...
        YYVPPGIPTNDTSNLER YYWNLSK
      rowData names(1): Sequence
      colnames(27): 12500amol_R1 12500amol_R2 ... 50amol_R2 50amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::rowData(D_norm_loess)
    Output
      DataFrame with 7944 rows and 1 column
                                              Sequence
                                           <character>
      AAAAQDEITGDGTTTVVCLVGELLR AAAAQDEITGDGTTTVVCLV..
      AAADALSDLEIK                        AAADALSDLEIK
      AAADALSDLEIKDSK                  AAADALSDLEIKDSK
      AAAEYEKGEYETAISTLNDAVEQGR AAAEYEKGEYETAISTLNDA..
      AAAPAQTTTDYK                        AAAPAQTTTDYK
      ...                                          ...
      YYTFNGPNYNENETIR                YYTFNGPNYNENETIR
      YYTITEVATR                            YYTITEVATR
      YYTLEEIQK                              YYTLEEIQK
      YYVPPGIPTNDTSNLER              YYVPPGIPTNDTSNLER
      YYWNLSK                                  YYWNLSK

---

    Code
      SummarizedExperiment::colData(D_norm_loess)
    Output
      DataFrame with 27 rows and 1 column
                         sample
                    <character>
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      125amol_R1     125amol_R1
      125amol_R2     125amol_R2
      ...                   ...
      500amol_R2     500amol_R2
      500amol_R3     500amol_R3
      50amol_R1       50amol_R1
      50amol_R2       50amol_R2
      50amol_R3       50amol_R3

---

    Code
      D_norm_lts
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(0):
      assays(1): logRatios
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(27): 12500amol_R1 12500amol_R2 ... 50amol_R2 50amol_R3
      colData names(1): sample

---

    Code
      SummarizedExperiment::rowData(D_norm_lts)
    Output
      DataFrame with 7944 rows and 1 column
                                              Sequence
                                           <character>
      AAAAQDEITGDGTTTVVCLVGELLR AAAAQDEITGDGTTTVVCLV..
      AAADALSDLEIK                        AAADALSDLEIK
      AAADALSDLEIKDSK                  AAADALSDLEIKDSK
      AAAEYEKGEYETAISTLNDAVEQGR AAAEYEKGEYETAISTLNDA..
      AAAPAQTTTDYK                        AAAPAQTTTDYK
      ...                                          ...
      YYTFNGPNYNENETIR                YYTFNGPNYNENETIR
      YYTITEVATR                            YYTITEVATR
      YYTLEEIQK                              YYTLEEIQK
      YYVPPGIPTNDTSNLER              YYVPPGIPTNDTSNLER
      YYWNLSK                                  YYWNLSK

---

    Code
      SummarizedExperiment::colData(D_norm_lts)
    Output
      DataFrame with 27 rows and 1 column
                         sample
                    <character>
      12500amol_R1 12500amol_R1
      12500amol_R2 12500amol_R2
      12500amol_R3 12500amol_R3
      125amol_R1     125amol_R1
      125amol_R2     125amol_R2
      ...                   ...
      500amol_R2     500amol_R2
      500amol_R3     500amol_R3
      50amol_R1       50amol_R1
      50amol_R2       50amol_R2
      50amol_R3       50amol_R3

---

    Code
      SummarizedExperiment::assays(D2)$logRatios
    Output
             logRatio_sample1_sample2 logRatio_sample1_sample3
      pep_1                        NA               1.09495817
      pep_2               -1.26147888               0.30166498
      pep_3                        NA               1.27132337
      pep_4               -0.52979568              -0.39790351
      pep_5                0.52218804               0.30487036
      pep_6                0.02686824              -0.18866120
      pep_7               -0.58551545              -1.64644581
      pep_8                0.11089975              -1.46836704
      pep_9               -1.44705619              -0.03532842
      pep_10               0.39241543              -1.00000000
             logRatio_sample2_sample3
      pep_1                 1.0949582
      pep_2                 1.5631439
      pep_3                 1.2713234
      pep_4                 0.1318922
      pep_5                -0.2173177
      pep_6                -0.2155294
      pep_7                -1.0609304
      pep_8                -1.5792668
      pep_9                 1.4117278
      pep_10               -1.3924154

---

    Code
      D2
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(1): imputed
      assays(2): logRatios maskImputation
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): logRatio_sample1_sample2 logRatio_sample1_sample3
        logRatio_sample2_sample3
      colData names(1): comparison

