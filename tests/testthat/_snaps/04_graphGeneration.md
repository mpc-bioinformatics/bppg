# test generateQuantGraphs

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]                
      [1,] "prot_4;prot_5" "pep_8;pep_9;pep_10"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
    Output
      [[1]]
      [1] NA NA
      
      [[2]]
      [1] 1.086 1.054 1.029
      

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
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
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
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]   
      [1,] "prot_3" "pep_5"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
    Output
      [[1]]
      [1] NA
      
      [[2]]
      [1] 0.964
      

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
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
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
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]         
      [1,] "prot_1" "pep_1"      
      [2,] "prot_1" "pep_3;pep_4"
      [3,] "prot_2" "pep_3;pep_4"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
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
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]         
      [1,] "prot_1" "pep_1;pep_2"
      [2,] "prot_1" "pep_3;pep_4"
      [3,] "prot_2" "pep_3;pep_4"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
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
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]         
      [1,] "prot_3"        "pep_5"      
      [2,] "prot_3"        "pep_7"      
      [3,] "prot_4;prot_5" "pep_7"      
      [4,] "prot_4;prot_5" "pep_8;pep_9"

---

    Code
      igraph::vertex_attr(graphs[[i]][[j]], "pep_ratio")
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
      

