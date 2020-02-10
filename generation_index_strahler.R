generation_index_strahler <- function(treeDF){
  ## This function counts the ordering in a tree based on the Horton-Strahler method.
  ## Strahler's stream order is a modification of Horton's stream order which fixes the ambiguity of Horton's ordering. In Strahler's ordering the main channel is not determined; instead the ordering is based on the hierarchy of tributaries. The ordering follows these rules:
  ## 1. if the node has no children, its Strahler order is 1.
  ## 2. if the node has one and only one tributary with Strahler greatest order i, and all other tributaries have order less than i, then the order remains i.
  ## 3. if the node has two or more tributaries with greatest order i, then the Strahler order of the node is i + 1.
  ##Strahler's stream ordering starts in initial links which assigns order one. It proceeds downstream. At every node it verifies that there are at least 2 equal tributaries with maximum order. If not, it continues with the highest order; if yes, it increases the node's order by 1 and continues downstream with the new order.

  generation <- mat.or.vec(length(treeDF$n),1)
  
  for (i in 1:length(treeDF$n)){
    if (is.na(treeDF$n[i])){
      generation[i] <- 1
    }
    else {generation[i] <- NA}
  }
  
  while(is.na(generation[1])){
    for (i in 1:length(treeDF$n)){
      sibs <- which(treeDF$parent == treeDF$nodeid[i])
      if (sum(is.na(generation[sibs])) == 0){
        for (j in sibs){
          parent <- which(treeDF$nodeid == treeDF$parent[j])
          if ((is.na(generation[parent])) & (sum(max(generation[sibs]) == generation[sibs]) >= 2)){
            generation[parent] <- max(generation[sibs]) + 1
          }
          else if (is.na(generation[parent])){
            generation[parent] <- max(generation[sibs])}
        }
      }
    }
  }
  ## works, but looks a bit weird for your highly asymmetric branches.  understimates in those cases.
  return(data.frame(treeDF,"generation" = generation))
}