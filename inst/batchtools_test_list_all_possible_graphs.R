library(batchtools)
library(data.table)

reg = makeExperimentRegistry(file.dir = "Registries/possible_bppgs",
                             packages = c("igraph", "rje", "rlist"))



problem_fun <- function(data, job, k) {
  return(NULL)
}

addProblem(name = "problem1", data = NULL, fun = problem_fun, reg = reg)

algo_fun <- function(data, instance, job, m) {

  ### Power set is set of possible nodes
  x <- rje::powerSet(1:m)
  x[[1]] <- NULL  ## delete empty set

  poss_node_comb <- t(rje::powerSetMat(length(x)))  ## alle möglichen Kombinationen von Knoten

  ### es muss mindestens log2(m+1) Knoten geben, sonst würden Proteine zusammengefasst
  poss_node_comb <- poss_node_comb[,colSums(poss_node_comb) >= log2(m+1)]

  G_list_tmp <- list()

  for (i in 1:ncol(poss_node_comb)) {
    #print(i)
    peptide_node_edges <- x[which(poss_node_comb[,i]==1)]

    prot_nodes <- unique(unlist(peptide_node_edges))
    if(length(prot_nodes) < m) next()

    M <- matrix(ncol = length(prot_nodes), nrow = length(peptide_node_edges))

    for(j in 1:length(peptide_node_edges)) {
      tmp <- prot_nodes %in% peptide_node_edges[[j]]
      M[j,] <- tmp
    }

    G <- igraph::graph_from_incidence_matrix(M)

    if (!igraph::is_connected(G)) next() ### nicht connected

    if(any(duplicated(M, MARGIN = 2))) next()  ### doppelte Spalten, Protein-Knoten würden gemergt

    ### isomorphic zu einem der schon bestehenden Graphen?

    if(length(G_list_tmp)) {
      test_iso <- sapply(G_list_tmp, function(x) {
        return(isomorphic_bipartite(x, G))
      })
      ## falls Graph nicht isomorph zu irgendeinem vorherigen ist, ergänze Liste
      if(!any(test_iso)) {
        G_list_tmp <- rlist::list.append(G_list_tmp, G)
      }
    } else {
      G_list_tmp <- rlist::list.append(G_list_tmp, G)
    }
  }
  return(G_list_tmp)
}

addAlgorithm(name = "generate_all_possible_bggps", fun = algo_fun)

ades <- list()
ades[[1]] <- data.table(m = 1:55)
names(ades) <- "generate_all_possible_bggps"

addExperiments(prob.designs = NULL, algo.designs = ades, repls = 1)


################################################################################

summarizeExperiments()

testJob(id = 2, external = TRUE)


reg$cluster.functions <- makeClusterFunctionsSSH(list(Worker$new("localhost", ncpus = 55)))


#chunks <- data.table(job.id = findJobs()$job.id, chunk = chunk(findJobs()$job.id, chunk.size = 10, shuffle = TRUE))

submitJobs(ids = findJobs()$job.id, resources = list(memory = 1024))

getStatus()










