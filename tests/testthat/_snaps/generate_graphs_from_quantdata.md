# test .generateQuantGraphs

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]          [,2]         
      [1,] "prot_1"      "pep_1;pep_2"
      [2,] "prot_1"      "pep_3;pep_4"
      [3,] "pep_3;pep_4" "prot_2"     

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]     [,2]   
      [1,] "prot_3" "pep_5"

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]                
      [1,] "prot_4;prot_5" "pep_10;pep_8;pep_9"

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]          [,2]         
      [1,] "prot_1"      "pep_1"      
      [2,] "prot_1"      "pep_3;pep_4"
      [3,] "pep_3;pep_4" "prot_2"     

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]                
      [1,] "prot_3"        "pep_6"             
      [2,] "prot_3"        "pep_7"             
      [3,] "pep_7"         "prot_4;prot_5"     
      [4,] "prot_4;prot_5" "pep_10;pep_8;pep_9"

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]          [,2]         
      [1,] "prot_1"      "pep_1;pep_2"
      [2,] "prot_1"      "pep_3;pep_4"
      [3,] "pep_3;pep_4" "prot_2"     

---

    Code
      igraph::as_edgelist(graphs[[i]][[j]])
    Output
           [,1]            [,2]           
      [1,] "prot_3"        "pep_5"        
      [2,] "prot_3"        "pep_7"        
      [3,] "pep_7"         "prot_4;prot_5"
      [4,] "prot_4;prot_5" "pep_8;pep_9"  

