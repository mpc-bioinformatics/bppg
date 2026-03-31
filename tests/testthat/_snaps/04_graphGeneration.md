# test generateQuantGraphs

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]                
      [1,] "prot_4;prot_5" "pep_8;pep_9;pep_10"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA NA
      
      [[2]]
      [1] 1.086 1.054 1.029
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA 1.029 1.086 1.054

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]         
      [1,] "prot_1" "pep_1;pep_2"
      [2,] "prot_1" "pep_3;pep_4"
      [3,] "prot_2" "pep_3;pep_4"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] NA
      
      [[3]]
      [1] 0.993 0.942
      
      [[4]]
      [1] 1.06 1.03
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA    NA 0.993 0.942 1.060 1.030

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]   
      [1,] "prot_3" "pep_5"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] 0.964
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA 0.964

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]       NA       NA 0.932000 0.900000 1.056333

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1]  TRUE  TRUE  TRUE  TRUE FALSE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] TRUE
      
      [[2]]
      [1] TRUE
      
      [[3]]
      [1] FALSE  TRUE
      
      [[4]]
      [1] TRUE
      
      [[5]]
      [1] FALSE FALSE FALSE
      

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]     NA     NA 0.9675 1.0450

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1] FALSE FALSE FALSE FALSE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] FALSE
      
      [[2]]
      [1] FALSE
      
      [[3]]
      [1] FALSE FALSE
      
      [[4]]
      [1] FALSE FALSE
      

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]                
      [1,] "prot_4;prot_5" "pep_7"             
      [2,] "prot_4;prot_5" "pep_8;pep_9;pep_10"
      [3,] "prot_3"        "pep_6"             
      [4,] "prot_3"        "pep_7"             

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA NA
      
      [[2]]
      [1] NA
      
      [[3]]
      [1] 1.086
      
      [[4]]
      [1] 0.9
      
      [[5]]
      [1] 0.953 0.955 1.004
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA    NA 1.004 1.086 0.900 0.953 0.955

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]         
      [1,] "prot_1" "pep_1"      
      [2,] "prot_1" "pep_3;pep_4"
      [3,] "prot_2" "pep_3;pep_4"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] NA
      
      [[3]]
      [1] 0.991
      
      [[4]]
      [1] 0.986 1.009
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA    NA 0.991 0.986 1.009

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]        NA        NA 0.9930000 0.9000000 0.9706667

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1] FALSE  TRUE  TRUE FALSE FALSE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] FALSE
      
      [[2]]
      [1] TRUE
      
      [[3]]
      [1]  TRUE FALSE
      
      [[4]]
      [1] FALSE
      
      [[5]]
      [1] FALSE FALSE FALSE
      

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]     NA     NA 0.9455 0.9975

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1]  TRUE FALSE  TRUE FALSE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] TRUE
      
      [[2]]
      [1] FALSE
      
      [[3]]
      [1] FALSE  TRUE
      
      [[4]]
      [1] FALSE FALSE
      

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]         
      [1,] "prot_1" "pep_1;pep_2"
      [2,] "prot_1" "pep_3;pep_4"
      [3,] "prot_2" "pep_3;pep_4"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] NA
      
      [[3]]
      [1] 0.945 0.982
      
      [[4]]
      [1] 1.023 0.943
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA    NA 0.945 0.982 1.023 0.943

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]         
      [1,] "prot_3"        "pep_5"      
      [2,] "prot_3"        "pep_7"      
      [3,] "prot_4;prot_5" "pep_7"      
      [4,] "prot_4;prot_5" "pep_8;pep_9"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] NA NA
      
      [[3]]
      [1] 1.033
      
      [[4]]
      [1] 0.92
      
      [[5]]
      [1] 1.022 1.006
      

---

    Code
      igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio")
    Output
      [1]    NA    NA 1.033 0.920 1.022 1.006

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]     NA     NA 0.9665 0.9200 0.9760

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1]  TRUE  TRUE  TRUE FALSE  TRUE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] TRUE
      
      [[2]]
      [1] TRUE
      
      [[3]]
      [1] FALSE  TRUE
      
      [[4]]
      [1] FALSE
      
      [[5]]
      [1] FALSE FALSE  TRUE
      

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "pep_ratio_mean")
    Output
      [1]     NA     NA 0.9635 0.9830

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "anyImputed")
    Output
      [1] FALSE FALSE FALSE FALSE

---

    Code
      igraph::vertex_attr(graphsImp[[i]][[j]], "imputed")
    Output
      [[1]]
      [1] FALSE
      
      [[2]]
      [1] FALSE
      
      [[3]]
      [1] FALSE FALSE
      
      [[4]]
      [1] FALSE FALSE
      

