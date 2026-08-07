# test .addUniquenessAttributes

    Code
      igraph::vertex_attr(G_new)
    Output
      $name
       [1] "sp|P00330|ADH1_YEAST"         "sp|P00331|ADH2_YEAST"        
       [3] "sp|P07246|ADH3_YEAST"         "sp|P38113|ADH5_YEAST"        
       [5] "ANELLINVK"                    "ANGTVVLVGLPAGAK"             
       [7] "ATDGGAHGVINVSVSEAAIEASTR"     "CCSDVFNQVVK"                 
       [9] "DIPVPEPKPNEILINVK"            "EALDFFAR"                    
      [11] "EALDFFSR"                     "EKDIVGAVLK"                  
      [13] "GVIFYESHGK"                   "IGDYAGIK"                    
      [15] "IQQGTDLAEVAPILCAGVTVYK"       "IVGLSELPK"                   
      [17] "LPLVGGHEGAGVVVGMGENVK"        "SANLMAGHWVAISGAAGGLGSLAVQYAK"
      [19] "SIGGEVFIDFTK"                 "SIPETQK"                     
      [21] "SISIVGSYVGNR"                 "VLGIDGGEGKEELFR"             
      [23] "VVGLSSLPEIYEK"                "VVGLSTLPEIYEK"               
      [25] "YSGVCHTDLHAWHGDWPLPVK"       
      
      $type
       [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
      [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
      [25] FALSE
      
      $pep_logRatio
       [1]          NA          NA          NA          NA -0.09213118  0.13908309
       [7] -0.63089964  0.18150780  0.10852732  0.02228366  0.07937947 -0.26744811
      [13]  0.21514238 -0.06718061  0.10167238  0.20413832 -0.21856356  0.02883078
      [19]  0.13946039  0.54554499 -0.06322847  0.22201082 -0.10181185 -0.03660714
      [25]  0.02723379
      
      $uniqueness
       [1]    NA    NA    NA    NA  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
      [13]  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE
      [25] FALSE
      
      $nr_unique_peptides
       [1]  9  2  4  0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
      
      $nr_shared_peptides
       [1]  6  5  1  1 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
      

# test .directBipartiteGraph

    Code
      igraph::as_edgelist(directed_graph)
    Output
           [,1]    [,2]    
      [1,] "pep_1" "prot_1"
      [2,] "pep_3" "prot_1"
      [3,] "pep_2" "prot_1"
      [4,] "pep_2" "prot_2"
      [5,] "pep_3" "prot_2"
      [6,] "pep_2" "prot_3"
      [7,] "pep_3" "prot_3"
      [8,] "pep_4" "prot_3"

---

    Code
      igraph::as_edgelist(directed_graph2)
    Output
           [,1]     [,2]   
      [1,] "prot_1" "pep_1"
      [2,] "prot_1" "pep_3"
      [3,] "prot_1" "pep_2"
      [4,] "prot_2" "pep_2"
      [5,] "prot_2" "pep_3"
      [6,] "prot_3" "pep_2"
      [7,] "prot_3" "pep_3"
      [8,] "prot_3" "pep_4"

