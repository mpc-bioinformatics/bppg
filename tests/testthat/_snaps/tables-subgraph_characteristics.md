# subgraph characteristics table

    Code
      res
    Output
        graph_ID nr_protein_nodes nr_peptide_nodes nr_unique_peptide_nodes
      1        1                7               13                       7
      2        2                1                1                       1
      3        3                2                2                       1
        nr_shared_peptide_nodes nr_edges nr_protein_accessions nr_peptide_sequences
      1                       6       22                     7                  476
      2                       0        1                     1                  204
      3                       1        3                     2                   47
        nr_prot_node_only_unique_pep nr_prot_node_unique_and_shared_pep
      1                            0                                  7
      2                            1                                  0
      3                            0                                  1
        nr_prot_node_only_shared_pep comparison
      1                            0          1
      2                            0          1
      3                            1          1

---

    Code
      res2
    Output
        graph_ID nr_protein_nodes nr_peptide_nodes nr_unique_peptide_nodes
      1        1                7               13                       7
      2        2                1                1                       1
      3        3                2                2                       1
      4        1                7               13                       7
      5        2                1                1                       1
      6        3                2                2                       1
        nr_shared_peptide_nodes nr_edges nr_protein_accessions nr_peptide_sequences
      1                       6       22                     7                  476
      2                       0        1                     1                  204
      3                       1        3                     2                   47
      4                       6       22                     7                  476
      5                       0        1                     1                  204
      6                       1        3                     2                   47
        nr_prot_node_only_unique_pep nr_prot_node_unique_and_shared_pep
      1                            0                                  7
      2                            1                                  0
      3                            0                                  1
      4                            0                                  7
      5                            1                                  0
      6                            0                                  1
        nr_prot_node_only_shared_pep comparison
      1                            0      comp1
      2                            0      comp1
      3                            1      comp1
      4                            0      comp2
      5                            0      comp2
      6                            1      comp2

