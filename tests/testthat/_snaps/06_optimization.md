# test .errorEquation

    Code
      

# test automated analysis

    Code
      res3
    Output
        Accession comparison graphID proteinNr error_optimal Ri_optimal Ci_optimal
      1        NA         NA      NA         1            NA         NA         NA
      2        NA         NA      NA         2            NA         NA         NA
           min_error error_constant  Ri       Ri_min    Ri_max Ci Ci_min Ci_max case
      1 8.326814e-18      partially  NA 8.323744e-08 0.5333333 NA    0.5    0.9    5
      2 8.326809e-18      partially 1.2           NA        NA NA    0.1    0.5    4

# test automated analysis with using results from other proteins

    Code
      res4
    Output
        Accession comparison graphID proteinNr error_optimal Ri_optimal Ci_optimal
      1        NA         NA      NA         1            NA         NA         NA
      2        NA         NA      NA         2            NA         NA         NA
           min_error error_constant  Ri       Ri_min    Ri_max Ci Ci_min Ci_max case
      1 8.326814e-18      partially  NA 8.323744e-08 0.5333333 NA    0.5    0.9    5
      2 8.326809e-18      partially 1.2           NA        NA NA    0.1    0.5    4

# test automated analysis with only one protein node

    Code
      res5
    Output
        Accession comparison graphID proteinNr error_optimal Ri_optimal Ci_optimal
      1        NA         NA      NA         1     0.2425386  0.8653497          1
        min_error error_constant        Ri Ri_min Ri_max Ci Ci_min Ci_max case
      1        NA             NA 0.8653497     NA     NA  1     NA     NA    6

# test automated analysis with constant error

    Code
      res6
    Output
        Accession comparison graphID proteinNr error_optimal Ri_optimal Ci_optimal
      1        NA         NA      NA         1            NA         NA         NA
      2        NA         NA      NA         2            NA         NA         NA
        min_error error_constant  Ri    Ri_min    Ri_max Ci Ci_min Ci_max case
      1     1e-05            yes 0.6        NA        NA NA    0.1    0.9    2
      2     1e-10            yes  NA 0.3523854 0.9752421 NA    0.1    0.9    1

