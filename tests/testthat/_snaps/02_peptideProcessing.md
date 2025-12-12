# test aggregateReplicates

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

# test calculatePeptideRatios

    Code
      SummarizedExperiment::assays(D1)$logRatios
    Output
             ratio_sample1_sample2 ratio_sample1_sample3 ratio_sample2_sample3
      pep_1                     NA                    NA           -0.01449257
      pep_2             0.23438359            0.30166498            0.06728139
      pep_3                     NA                    NA           -0.13742681
      pep_4            -0.52979568           -0.39790351            0.13189218
      pep_5             0.52218804            0.30487036           -0.21731768
      pep_6             0.02686824           -0.18866120           -0.21552944
      pep_7            -0.58551545                    NA                    NA
      pep_8             0.11089975           -0.26790101           -0.37880076
      pep_9            -0.11350943           -0.03532842            0.07818101
      pep_10            0.39241543            0.44831953            0.05590409

---

    Code
      D1
    Output
      class: SummarizedExperiment 
      dim: 10 3 
      metadata(0):
      assays(1): logRatios
      rownames(10): pep_1 pep_2 ... pep_9 pep_10
      rowData names(1): Sequence
      colnames(3): ratio_sample1_sample2 ratio_sample1_sample3
        ratio_sample2_sample3
      colData names(1): comparison

