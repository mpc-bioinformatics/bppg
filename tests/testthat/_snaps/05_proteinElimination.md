# test proteinElimination

    Code
      res
    Output
      $min_error_ref
      [1] 0.02231625
      
      $protsOriginIDs
      [1] "prot_1" "prot_2" "prot_3" "prot_4" "prot_5"
      
      $resDF
               comb n_proteins   error current_best
      1  1,5,6,8,10          5 0.02232         TRUE
      2    5,6,8,10          4 0.02232         TRUE
      3      6,8,10          3      NA        FALSE
      4      5,8,10          3 0.02232         TRUE
      5        8,10          2      NA        FALSE
      6        5,10          2 0.05975        FALSE
      7         5,8          2      NA        FALSE
      8      5,6,10          3 0.05975        FALSE
      9       5,6,8          3 0.02232         TRUE
      10        6,8          2      NA        FALSE
      11        5,6          2      NA        FALSE
      12   1,6,8,10          4 0.02232        FALSE
      13     1,8,10          3 0.02232         TRUE
      14       1,10          2 0.02232         TRUE
      15         10          1      NA        FALSE
      16          1          1      NA        FALSE
      17        1,8          2 0.04905        FALSE
      18     1,6,10          3 0.02232        FALSE
      19       6,10          2      NA        FALSE
      20        1,6          2      NA        FALSE
      21      1,6,8          3 0.03459        FALSE
      22   1,5,8,10          4 0.02232        FALSE
      23     1,5,10          3 0.02232        FALSE
      24        1,5          2      NA        FALSE
      25      1,5,8          3 0.02232        FALSE
      26   1,5,6,10          4 0.02232        FALSE
      27      1,5,6          3      NA        FALSE
      28    1,5,6,8          4 0.00788        FALSE
      
      $res_best
      $res_best$G
      $res_best$G[[1]]
           [,1]     [,2]    
      [1,] "prot_1" "pep_1" 
      [2,] "prot_1" "pep_2" 
      [3,] "prot_1" "pep_3" 
      [4,] "pep_3"  "prot_5"
      [5,] "pep_4"  "prot_5"
      [6,] "pep_5"  "prot_5"
      
      
      $res_best$comb
      [1] "1,10"
      
      $res_best$n_comb
      [1] 2
      
      $res_best$error
      [1] 0.02231625
      
      

