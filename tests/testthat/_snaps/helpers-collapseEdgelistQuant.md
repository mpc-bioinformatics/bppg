# test .collapseEdgelistQuant

    Code
      collapsed_edgelist
    Output
                     protein     peptide pep_ratio
      1               prot_1 pep_1;pep_2  0.9;1.02
      2               prot_2 pep_1;pep_2  0.9;1.02
      5               prot_2       pep_3      0.96
      6 prot_3;prot_4;prot_5       pep_3      0.96
      9 prot_3;prot_4;prot_5 pep_4;pep_5 0.96;1.06

---

    Code
      collapsed_edgelist2
    Output
                      protein peptide pep_ratio
      1                prot_1   pep_1      1.02
      2                prot_2   pep_1      1.02
      3                prot_1   pep_2      0.90
      4                prot_2   pep_2      0.90
      5                prot_2   pep_3      0.96
      6  prot_3;prot_4;prot_5   pep_3      0.96
      9  prot_3;prot_4;prot_5   pep_4      0.96
      12 prot_3;prot_4;prot_5   pep_5      1.06

