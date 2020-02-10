### This script explores consequences of using Horton Strahler method for labeling generations and potential implications it has for length scaling.  It is being tested out on fully referenced artery/vein data from Charbonnier et al., reconstructed using Angicart.

library("rgl")
library("mgcv")
library("nat")
library("geometry")

source("/Users/alexwork/odrive/Dropbox/Research/UA/GRRRangerism_online/asymmetric_modelling/R_functions/generation_index_strahler.R")

setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/lung_cancer_project/carve14/av_reconstruction/full_referenced_output/2064.1140250077/")

arteries <- read.csv(file = "arteries/2064.1140250077_arteries_resorted.csv", header = TRUE)
veins <- read.csv(file = "veins/2064.1140250077_veins_resorted.csv", header = TRUE)

# Grab single large trees.
art_tree_big <- arteries[which(arteries$treeid == 1),]
vein_tree_big <- veins[which(veins$treeid == 4),] # using tree 4 as it is largest with best identified root.

nopen3d()
plot3d(arteries[,1:3], col = "blue")
plot3d(veins[,1:3], col = "red", add = TRUE)

# Graph to examine trees
nopen3d()
plot3d(art_tree_big[,1:3])
plot3d(art_tree_big[which(is.na(art_tree_big$parent)),1:3], size = 5, col = "green", add = TRUE)

nopen3d()
plot3d(vein_tree_big[,1:3])
plot3d(vein_tree_big[which(is.na(vein_tree_big$parent)), 1:3], size = 5, col = "green", add = TRUE)


# Add Horton-Strahler labeling.  Need to first make unique matrix sans xyz coords, then paste with replication for xyz coords.

unique_art_tree_big <- unique(art_tree_big[, c("nodeid", "parent", "n")])
unique_art_tree_big <- generation_index_strahler(unique_art_tree_big)
nonunique_art_tre_big_hsgen <- mat.or.vec(nrow(art_tree_big), 1)
for(i in 1:length(unique_art_tree_big$nodeid)){
  rep_rows <- which(art_tree_big$nodeid == unique_art_tree_big$nodeid[i])
  nonunique_art_tre_big_hsgen[rep_rows] <- unique_art_tree_big$generation[i]
}
art_tree_big <- data.frame(art_tree_big, "generation" = nonunique_art_tre_big_hsgen)

unique_vein_tree_big <- unique(vein_tree_big[, c("nodeid", "parent", "n")])
unique_vein_tree_big <- generation_index_strahler(unique_vein_tree_big)
nonunique_vein_tre_big_hsgen <- mat.or.vec(nrow(vein_tree_big), 1)
for(i in 1:length(unique_vein_tree_big$nodeid)){
  rep_rows <- which(vein_tree_big$nodeid == unique_vein_tree_big$nodeid[i])
  nonunique_vein_tre_big_hsgen[rep_rows] <- unique_vein_tree_big$generation[i]
}
vein_tree_big <- data.frame(vein_tree_big, "generation" = nonunique_vein_tre_big_hsgen)


## Graph the large arterial and venous networks using original node labeling as coloring and the Horton-Strahler labeling as coloring for comparison.


nopen3d()
plot3d(art_tree_big[,1:3], col = art_tree_big$generation)
text3d(art_tree_big[,1:3], texts = art_tree_big$generation, add = TRUE)

nopen3d()
plot3d(vein_tree_big[,1:3], col = vein_tree_big$nodeid+1)
text3d(vein_tree_big[,1:3], texts = vein_tree_big$generation, add = TRUE)


art_chain_list <- list()

## The following loop should correctly identify which vessels comprise distinct vessel-chains of a given generation, labeled by name.  output is a list of lists, where the primary list indexing corresponds to the horton-strahler indexing, the secondary listing corresponds to distint vessel-chains for a given horton strahler value, then within these lists are vectors of vessel names in a given chain.

for(i in sort(unique(art_tree_big$generation), decreasing = TRUE)){
  subtree <- unique(art_tree_big[which(art_tree_big$generation == i),c("nodeid", "parent", "generation")])
  art_chain_list[[i]] <- list(c(subtree$nodeid[1]))
  labeled_vessels <- c(subtree$nodeid[1])
  for(j in 2:nrow(subtree)){
    for(k in 1:length(art_chain_list[[i]])){
      if(subtree$parent[j] %in% art_chain_list[[i]][[k]]){
        art_chain_list[[i]][[k]] <- c(art_chain_list[[i]][[k]], subtree$nodeid[j])
        labeled_vessels <- c(labeled_vessels, subtree$nodeid[j])
      }else if(!(subtree$parent[j] %in% art_chain_list[[i]][[k]]) & !((subtree$parent[j] %in% labeled_vessels) | (subtree$nodeid[j] %in% labeled_vessels))){
        art_chain_list[[i]][[length(art_chain_list[[i]])+1]] <- c(subtree$nodeid[j])
        labeled_vessels <- c(labeled_vessels, subtree$nodeid[j])
      }
    }
  }
}

# here we plot vessel chain by vessel chain.  manual inspection so far looks good.  Consider using art_chain_list to find volumes/lengths/etc. and add all together and compare to sum from treeDF.

colors <- list(rep("black", length(art_chain_list[[1]])), rep("red", length(art_chain_list[[2]])), rep("green", length(art_chain_list[[3]])))

nopen3d()
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[4]][[1]]), 1:3], col = "blue")
for(i in 1:4){
  for(j in 1:length(art_chain_list[[i]])){
    plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[i]][[j]]), 1:3], col = colors[[i]][j], add = TRUE)
  }
}
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[1]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[2]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[3]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[4]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[5]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[3]][[6]]), 1:3], col = "green", add = TRUE)
plot3d(art_tree_big[which(art_tree_big$nodeid %in% art_chain_list[[2]][[1]]), 1:3], col = "red", add = TRUE)


art_chain_list

## Now to relabel vessels based on their Horton-Strahler labeling.  Each vector in each list can be given a unique, and repeated id,  call it hs_nodeid.  The challenging part will be identifying and relabeling parents, call them hs_parent.

# split the original nodeid to extract the unique treeid.  Or, use the treeid value from the treeid column!
count <- 1
treeid <- strsplit(x = as.character(art_tree_big$nodeid[1]), split = "\\.")[[1]][2] # note this is a character to be pasted


hs_art_chain_list <- list()
for(i in length(art_chain_list):1){
  hs_art_chain_list[[i]] <- list()
  for(j in 1:length(art_chain_list[[i]])){
    hs_art_chain_list[[i]][[j]] <- as.numeric(paste(count, treeid, sep = "."))
    count <- count + 1
  }
}

# Now we can map from old vessel ids to new vessel ids, and doing so using rows in original data frame.  This means hs_nodeids now serve a key for mapping.  Should be portable to parents as well as measured values, except the lengths need to be added, ?and the radii averaged?
hs_nodeid <- mat.or.vec(nrow(art_tree_big), 1)
hs_temp_parent <- mat.or.vec(nrow(art_tree_big), 1)
hs_parent <- mat.or.vec(nrow(art_tree_big), 1)

for(i in 1:length(art_chain_list)){
  # i <- 3
  for(j in 1:length(art_chain_list[[i]])){
    # j <- 1
    rows <- which(art_tree_big$nodeid %in% art_chain_list[[i]][[j]])
    hs_nodeid[rows] <- hs_art_chain_list[[i]][[j]]
    hs_temp_parent[rows] <- rep(unique(art_tree_big$parent[rows][which(!(art_tree_big$parent[rows] %in% art_tree_big$nodeid[rows]))]), length(rows))
    # Find hs_nodeid corresponding to this newly found parent value (from old labeling), then reassign.  Considering using this current assignment for a temporary parent then jumping to the key.  Do if not an NA value too for the root.
  }
}

for(i in 1:length(hs_temp_parent)){
  if(!is.na(hs_temp_parent[i])){
    # i <- 10
    hs_parent[i] <- unique(hs_nodeid[which(art_tree_big$nodeid == hs_temp_parent[i])])
  }else
    hs_parent[i] <- NA
}

## Consider keeping radii so we have local measure since tapering occurs.  For lengths, sample code for calculating total length is below.  Use this to assign new repeated length value an hs_length column for all instances where hs_nodeid == 1.001, then loop through all hs_nodeids.
hs_length <- mat.or.vec(nrow(art_tree_big),1)

for(i in 1:length(unique(hs_nodeid))){
  vessel <- unique(hs_nodeid)[i]
  rows <- which(hs_nodeid == vessel)
  new_length <- sum(unique(art_tree_big$length[rows]))
  hs_length[rows] <- new_length
}

# Merge to main art_tree_big_object
art_tree_big <- data.frame(art_tree_big, "hs_nodeid" = hs_nodeid, "hs_parent" = hs_parent, "hs_length" = hs_length)

# Make seperate unique art_tree_big for HS quantities
hs_unique_art_tree_big <- unique(art_tree_big[ , c("generation", "hs_nodeid", "hs_parent", "hs_length")])

## Quickly examine length ratio using this approach.

source("/Users/alexwork/odrive/Dropbox/Research/UA/GRRRangerism_online/asymmetric_modelling/R_functions/get.ratio.R")

hist(get.ratio(treeDF = hs_unique_art_tree_big, treelist = hs_unique_art_tree_big$hs_length, is.hs = TRUE), xlim = c(0, 5), breaks = 100)



nopen3d()
plot3d(art_tree_big[,1:3], col = art_tree_big$generation)
text3d(art_tree_big[,1:3], texts = art_tree_big$hs_nodeid, add = TRUE)



# Generalize over all subtrees now by writing as one function

horton_strahler_relabel <- function(treeDF){
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
  treeid <- strsplit(x = as.character(treeDF$nodeid[1]), split = "\\.")[[1]][2] # note this is a character to be pasted
  
  hs_chain_list <- list()
  for(i in length(chain_list):1){
    hs_chain_list[[i]] <- list()
    for(j in 1:length(chain_list[[i]])){
      hs_chain_list[[i]][[j]] <- as.numeric(paste(count, treeid, sep = "."))
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
  
  # Merge to main art_tree_big_object
  treeDF <- data.frame(treeDF, "hs_nodeid" = hs_nodeid, "hs_parent" = hs_parent, "hs_length" = hs_length)
  
  return(treeDF)
}

# Now we can run over all trees in the arteries and veins data frame.

# Need to create a new dataframe to save output to as have to process each sub tree individually.

# hs_arteries <- arteries[0,]
# hs_arteries <- data.frame(hs_arteries, "generation" = c(NULL), "hs_nodeid" = c(NULL), "hs_parent" = c(NULL), "hs_length" = c(NULL))

# Identify unique tree id values to loop through
treeids <- unique(veins$treeid)
# Initialize the final data frame by running horton_strahler_relabel on first tree id from treeids vector.
hs_veins <- horton_strahler_relabel(treeDF = veins[which(veins$treeid == treeids[1]), ])
# Loop through remaining treeids.
for(i in treeids[-1]){
  hs_veins <- rbind(hs_veins, horton_strahler_relabel(treeDF = veins[which(veins$treeid == i), ]))
}

## Quickly examine length ratio using this approach.
source("/Users/alexwork/odrive/Dropbox/Research/UA/GRRRangerism_online/asymmetric_modelling/R_functions/get.ratio.R")

hs_veins_unique <- unique(hs_veins[ , c("generation", "gen_diff", "hs_nodeid", "hs_parent", "hs_length")])
veins_unique <- unique(veins[ , c("nodeid", "parent", "length")])

gamma <- get.ratio(treeDF = veins_unique, treelist = veins_unique$length)

hs_gamma <- get.ratio(treeDF = hs_veins_unique, treelist = hs_veins_unique$hs_length, is.hs = TRUE)
hs_veins_unique <- data.frame(hs_veins_unique, "hs_gamma" = hs_gamma)

hist(hs_gamma, breaks = 500, xlim = c(0, 20))
hist(gamma, breaks = 500, xlim = c(0, 20))

median(gamma, na.rm = T)
median(hs_gamma, na.rm = T)

hist(hs_gamma[which(hs_veins_unique$gen_diff == 1)], breaks = 500, xlim = c(0,5))
median(hs_gamma[which(hs_veins_unique$gen_diff == 1)], na.rm = T)

hist(hs_gamma[which(hs_veins_unique$gen_diff == 2)], breaks = 500, xlim = c(0,5))
median(hs_gamma[which(hs_veins$gen_diff == 2)], na.rm = T)

hist(hs_gamma[which(hs_veins_unique$gen_diff == 3)], breaks = 500, xlim = c(0,5))
median(hs_gamma[which(hs_veins_unique$gen_diff == 3)], na.rm = T)

## Loop over all arteries and veins that you have.
## Try for the plant data that you have access to.
## Also, bin by different integer differences between branches (i + 1, i + 2, i + 3, ...)

