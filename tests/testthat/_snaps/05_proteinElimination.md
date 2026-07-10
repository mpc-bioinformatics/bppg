# test proteinElimination

    Code
      res
    Output
      $min_error_ref
      [1] 1.005431
      
      $protsOriginIDs
      [1] "sp|P00330|ADH1_YEAST" "sp|P00331|ADH2_YEAST" "sp|P07246|ADH3_YEAST"
      [4] "sp|P38113|ADH5_YEAST"
      
      $resDF
           comb n_proteins  error current_best
      1 1,2,3,4          4 1.0054         TRUE
      2   2,3,4          3     NA        FALSE
      3   1,3,4          3     NA        FALSE
      4   1,2,4          3     NA        FALSE
      5   1,2,3          3 1.0056         TRUE
      6     2,3          2     NA        FALSE
      7     1,3          2     NA        FALSE
      8     1,2          2     NA        FALSE
      
      $res_best
      $res_best$G
      $res_best$G[[1]]
            [,1]                   [,2]                          
       [1,] "sp|P00330|ADH1_YEAST" "ANELLINVK"                   
       [2,] "sp|P00331|ADH2_YEAST" "ANGTVVLVGLPAGAK"             
       [3,] "sp|P00330|ADH1_YEAST" "ATDGGAHGVINVSVSEAAIEASTR"    
       [4,] "sp|P00330|ADH1_YEAST" "CCSDVFNQVVK"                 
       [5,] "sp|P07246|ADH3_YEAST" "DIPVPEPKPNEILINVK"           
       [6,] "sp|P00330|ADH1_YEAST" "EALDFFAR"                    
       [7,] "sp|P00331|ADH2_YEAST" "EALDFFAR"                    
       [8,] "sp|P07246|ADH3_YEAST" "EALDFFSR"                    
       [9,] "sp|P00330|ADH1_YEAST" "EKDIVGAVLK"                  
      [10,] "sp|P00330|ADH1_YEAST" "GVIFYESHGK"                  
      [11,] "sp|P00330|ADH1_YEAST" "IGDYAGIK"                    
      [12,] "sp|P00331|ADH2_YEAST" "IGDYAGIK"                    
      [13,] "sp|P07246|ADH3_YEAST" "IQQGTDLAEVAPILCAGVTVYK"      
      [14,] "sp|P07246|ADH3_YEAST" "IVGLSELPK"                   
      [15,] "sp|P00330|ADH1_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [16,] "sp|P00331|ADH2_YEAST" "LPLVGGHEGAGVVVGMGENVK"       
      [17,] "sp|P00330|ADH1_YEAST" "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [18,] "sp|P00330|ADH1_YEAST" "SIGGEVFIDFTK"                
      [19,] "sp|P00330|ADH1_YEAST" "SIPETQK"                     
      [20,] "sp|P00331|ADH2_YEAST" "SIPETQK"                     
      [21,] "sp|P00330|ADH1_YEAST" "SISIVGSYVGNR"                
      [22,] "sp|P00331|ADH2_YEAST" "SISIVGSYVGNR"                
      [23,] "sp|P00330|ADH1_YEAST" "VLGIDGGEGKEELFR"             
      [24,] "sp|P00331|ADH2_YEAST" "VVGLSSLPEIYEK"               
      [25,] "sp|P00330|ADH1_YEAST" "VVGLSTLPEIYEK"               
      [26,] "sp|P00330|ADH1_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      [27,] "sp|P07246|ADH3_YEAST" "YSGVCHTDLHAWHGDWPLPVK"       
      
      
      $res_best$comb
      [1] "1,2,3"
      
      $res_best$n_comb
      [1] 3
      
      $res_best$error
      [1] 1.00567
      
      

