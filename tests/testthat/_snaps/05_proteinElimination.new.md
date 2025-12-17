# test proteinElimination

    Code
      res
    Output
      $min_error_ref
      [1] 0.01
      
      $min_error_current
      [1] 0.01
      
      $protein_nodes_list
      [1] "prot_1" "prot_2" "prot_3" "prot_4" "prot_5"
      
      $combination_list
       [1] "1,5,6,8,10" "5,6,8,10"   "6,8,10"     "5,8,10"     "8,10"      
       [6] "5,10"       "5,8"        "5,6,10"     "5,6,8"      "6,8"       
      [11] "5,6"        "1,6,8,10"   "1,5,8,10"   "1,5,6,10"   "1,5,6,8"   
      
      $error_list
       [1] 0.010 0.010    NA 0.010    NA 0.025    NA 0.025 0.010    NA    NA    NA
      [13]    NA    NA    NA
      
      $comb_current
      [1] "5,6,8"
      
      $G_current
           [,1]     [,2]    
      [1,] "pep_1"  "prot_2"
      [2,] "pep_2"  "prot_2"
      [3,] "pep_2"  "prot_3"
      [4,] "pep_3"  "prot_3"
      [5,] "prot_3" "pep_4" 
      [6,] "pep_4"  "prot_4"
      [7,] "prot_4" "pep_5" 
      
      $n_comb_current
      [1] 3
      

