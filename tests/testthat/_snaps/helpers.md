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
       [1]          NA          NA          NA          NA -0.09208193  0.15016353
       [7] -0.63333305  0.19277551  0.09378775  0.02149638  0.07320356 -0.25189929
      [13]  0.21497370 -0.06578712  0.09921028  0.18956323 -0.21712869  0.03261818
      [19]  0.14143327  0.53738274 -0.06636144  0.21885569 -0.08659069 -0.03884356
      [25]  0.02805596
      
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

