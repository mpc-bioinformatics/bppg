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

