library(batchtools)
library(data.table)


install_github("d4ndo/binaryLogic")

reg = makeExperimentRegistry(file.dir = "Registries/possible_bppgs",
                             source = c("promotion_project/R scripts/helpers-isomorphisms.R"),
                             packages = c("igraph", "rje", "rlist", "binaryLogic"))



problem_fun <- function(data, job, k) {
  return(NULL)
}

addProblem(name = "problem1", data = NULL, fun = problem_fun, reg = reg)

algo_fun <- function(data, instance, job, m) {

  ### Power set is set of possible nodes
  x <- rje::powerSet(1:m)
  x[[1]] <- NULL  ## delete empty set
  #x

  nr_poss_pep_nodes <- length(x)

  nr_poss_node_comb <- 2^nr_poss_pep_nodes-1 # wieder Potenzmenge aller möglichen Knoten (außer leere Menge)
  ### leere Menge wird dann so codiert

  G_list_tmp <- list()

  for (i in 1:nr_poss_node_comb) {
    print(i)

    ind <- binaryLogic::as.binary(i, n = nr_poss_pep_nodes)

    peptide_node_edges <- x[ind]

    ### TODO: Abfragen, dass alle Proteinknoten erwischt werden
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
ades[[1]] <- data.table(m = 1:6)
names(ades) <- "generate_all_possible_bggps"

addExperiments(prob.designs = NULL, algo.designs = ades, repls = 1)


################################################################################

summarizeExperiments()

testJob(id = 2, external = TRUE)


reg$cluster.functions <- makeClusterFunctionsSSH(list(Worker$new("localhost", ncpus = 55)))


#chunks <- data.table(job.id = findJobs()$job.id, chunk = chunk(findJobs()$job.id, chunk.size = 10, shuffle = TRUE))

submitJobs(ids = findJobs()$job.id, resources = list(memory = 1024))

getStatus()










