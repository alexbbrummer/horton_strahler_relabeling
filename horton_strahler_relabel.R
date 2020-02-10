horton_strahler_relabel <- function(treeDF, treeid = NULL){
  ## This function relabels a tree with horton strahler labeling, and merges branches of similar generation into one single branch, thus lengthening the branch. Raddi are maintained as is to track tapering during side-branching.  The original tree data frame is returned with the newly labeled horton strahler generations, node ids, parent ids, and lengths appended.
  source("/Users/alexwork/odrive/Dropbox/Research/UA/GRRRangerism_online/asymmetric_modelling/R_functions/generation_index_strahler.R")
  
  ## First we add the Horton-Strahler labelling.  Need to first make a unique matrix without the xyz coordinates, then paste with replication the generation labeling back onto the original data frame.  Note that if we are working with trees that don't have xyz coordinates then we may need to provide a use-case parameter in the function to toggle this routine.
  
  unique_treeDF <- unique(treeDF[, c("nodeid", "parent", "n")])
  unique_treeDF <- generation_index_strahler(unique_treeDF)
  nonunique_treeDF_hsgen <- mat.or.vec(nrow(treeDF), 1)
  for(i in 1:length(unique_treeDF$nodeid)){
    rep_rows <- which(treeDF$nodeid == unique_treeDF$nodeid[i])
    nonunique_treeDF_hsgen[rep_rows] <- unique_treeDF$generation[i]
  }
  treeDF <- data.frame(treeDF, "generation" = nonunique_treeDF_hsgen)
  
  ## The following loop should correctly identify which vessels comprise distinct vessel-chains of a given generation, labeled by name.  Output is a list of lists, where the primary list indexing corresponds to the horton-strahler indexing, the secondary listing corresponds to distint vessel-chains for a given horton strahler value, then within these lists are vectors of vessel names in a given chain.
  
  chain_list <- list()
  
  for(i in sort(unique(treeDF$generation), decreasing = TRUE)){
    # Loop first by generation, and subset by all vessels of common generation.
    subtree <- unique(treeDF[which(treeDF$generation == i),c("nodeid", "parent", "generation")])
    # Initialize chain list
    chain_list[[i]] <- list(c(subtree$nodeid[1]))
    # Populate chain list with first vessel entry.
    labeled_vessels <- c(subtree$nodeid[1])
    for(j in 2:nrow(subtree)){
      # Initialize secondary chain list
      for(k in 1:length(chain_list[[i]])){
        # Check if parent of vessel j exists in list.  If so, add vessel j to same list, and add to vector of labeled vessels.
        if(subtree$parent[j] %in% chain_list[[i]][[k]]){
          chain_list[[i]][[k]] <- c(chain_list[[i]][[k]], subtree$nodeid[j])
          labeled_vessels <- c(labeled_vessels, subtree$nodeid[j])
          # If parent of vessel j is not in list, and it doesn't exist in a different list as indicated by existing in the vector of labeled vessels nor vessel j is not in the vector of labeled vessels, then add vessel j to both chain list and labeled vessels.
        }else if(!(subtree$parent[j] %in% chain_list[[i]][[k]]) & !((subtree$parent[j] %in% labeled_vessels) | (subtree$nodeid[j] %in% labeled_vessels))){
          chain_list[[i]][[length(chain_list[[i]])+1]] <- c(subtree$nodeid[j])
          labeled_vessels <- c(labeled_vessels, subtree$nodeid[j])
        }
      }
    }
  }
  
  
  ## Now to relabel vessels based on their Horton-Strahler labeling.  Each vector in each list can be given a unique, and repeated, id, call it hs_nodeid.  The hs parent vector will be labeled as hs_parent.
  
  # split the original nodeid to extract the unique treeid.  Or, use the treeid value from the treeid column.  This is particular for the angicart++ output tree files.
  count <- 1
  if(is.null(treeid)){
    treeid <- strsplit(x = as.character(treeDF$nodeid[1]), split = "\\.")[[1]][2] # note this is a character to be pasted
  }
  
  
  hs_chain_list <- list()
  for(i in length(chain_list):1){
    hs_chain_list[[i]] <- list()
    for(j in 1:length(chain_list[[i]])){
      if(is.null(treeid)){
        hs_chain_list[[i]][[j]] <- as.numeric(paste(count, treeid, sep = "."))
      }else{
        hs_chain_list[[i]][[j]] <- paste(count, treeid, sep = ".")
      }
      count <- count + 1
    }
  }
  
  # Now we can map from old vessel ids to new vessel ids, and doing so using rows in original data frame.  This means hs_nodeids now serve a key for mapping.
  
  hs_nodeid <- mat.or.vec(nrow(treeDF), 1)
  hs_temp_parent <- mat.or.vec(nrow(treeDF), 1)
  hs_parent <- mat.or.vec(nrow(treeDF), 1)
  
  for(i in 1:length(chain_list)){
    for(j in 1:length(chain_list[[i]])){
      rows <- which(treeDF$nodeid %in% chain_list[[i]][[j]])
      hs_nodeid[rows] <- hs_chain_list[[i]][[j]]
      # Here we identify which vessels are associated with a given list-list entry, and what their parent ids are.  Then, the one parent id not included in the vessel id list gives us a global parent for the whole chain.  This parent id is used in next for loop to identify new parent according to chain-list key. 
      hs_temp_parent[rows] <- rep(unique(treeDF$parent[rows][which(!(treeDF$parent[rows] %in% treeDF$nodeid[rows]))]), length(rows))
    }
  }
  
  # Find hs_nodeid corresponding to this newly found parent value (from old labeling), then reassign.  Considering using this current assignment for a temporary parent then jumping to the key.  If !is.na is used to preserve the root-parent label.
  
  for(i in 1:length(hs_temp_parent)){
    if(!is.na(hs_temp_parent[i])){
      hs_parent[i] <- unique(hs_nodeid[which(treeDF$nodeid == hs_temp_parent[i])])
    }else
      hs_parent[i] <- NA
  }
  
  ## Here we calculate the new total lengths
  
  hs_length <- mat.or.vec(nrow(treeDF),1)
  
  for(i in 1:length(unique(hs_nodeid))){
    vessel <- unique(hs_nodeid)[i]
    rows <- which(hs_nodeid == vessel)
    new_length <- sum(unique(treeDF$length[rows]))
    hs_length[rows] <- new_length
  }
  
  ## Here we calculate the difference in generation between parent and child branches.
  
  gen_diff <- mat.or.vec(nrow(treeDF), 1)
  gen_diff[] <- NA
  for(i in 1:length(gen_diff)){
    if(!is.na(hs_parent[i])){
      gen_diff[i] <- unique(treeDF$generation[which(hs_nodeid == hs_parent[i])]) - treeDF$generation[i]
    }
  }
  
  # Merge to main art_tree_big_object
  treeDF <- data.frame(treeDF, "gen_diff" = gen_diff, "hs_nodeid" = as.numeric(as.character(hs_nodeid)), "hs_parent" = as.numeric(as.character(hs_parent)), "hs_length" = hs_length)
  
  return(treeDF)
}