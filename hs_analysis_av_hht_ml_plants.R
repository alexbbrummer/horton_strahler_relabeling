## This script examines how the Horton Strahler relabeling changes length scale factor estimates.  We will explore changes in the humna lung arteries and veins data, as well as the human head and torso, mouse lung, and plant data from earlier papers.

source("/Users/alexwork/odrive/Dropbox/Research/UA/GRRRangerism_online/asymmetric_modelling/R_functions/get.ratio.R")
source("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/lung_cancer_project/carve14/horton_strahler_relabel.R")

#### Reading in and recalculating labeling and scale factors for AV data. ####
setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/lung_cancer_project/carve14/av_reconstruction/full_referenced_output/")

# Define vectors of filenames and pathnames
filenamepath_vec_arteries_resorted <- c("1396.1132404220/arteries/1396.1132404220_arteries_resorted.csv", "2064.1140250077/arteries/2064.1140250077_arteries_resorted.csv", "2604.1126357612/arteries/2604.1126357612_arteries_resorted.csv", "3300.1120308766/arteries/3300.1120308766_arteries_resorted.csv", "3356.1097932728/arteries/3356.1097932728_arteries_resorted.csv", "3476.1102167727/arteries/3476.1102167727_arteries_resorted.csv", "3840.1105186168/arteries/3840.1105186168_arteries_resorted.csv", "4332.1132400990/arteries/4332.1132400990_arteries_resorted.csv", "4584.1131197859/arteries/4584.1131197859_arteries_resorted.csv", "4840.1113050391/arteries/4840.1113050391_arteries_resorted.csv")

filenamepath_vec_veins_resorted <- c("1396.1132404220/veins/1396.1132404220_veins_resorted.csv", "2064.1140250077/veins/2064.1140250077_veins_resorted.csv", "2604.1126357612/veins/2604.1126357612_veins_resorted.csv", "3300.1120308766/veins/3300.1120308766_veins_resorted.csv", "3356.1097932728/veins/3356.1097932728_veins_resorted.csv", "3476.1102167727/veins/3476.1102167727_veins_resorted.csv", "3840.1105186168/veins/3840.1105186168_veins_resorted.csv", "4332.1132400990/veins/4332.1132400990_veins_resorted.csv", "4584.1131197859/veins/4584.1131197859_veins_resorted.csv", "4840.1113050391/veins/4840.1113050391_veins_resorted.csv")

filenamepath_vec_arteries <- c("1396.1132404220/arteries/1396.1132404220_arteries_withRoots.tsv", "2064.1140250077/arteries/2064.1140250077_arteries_withRoots.tsv", "2604.1126357612/arteries/2604.1126357612_arteries_withRoots.tsv", "3300.1120308766/arteries/3300.1120308766_arteries_withRoots.tsv", "3356.1097932728/arteries/3356.1097932728_arteries_withRoots.tsv", "3476.1102167727/arteries/3476.1102167727_arteries_withRoots.tsv", "3840.1105186168/arteries/3840.1105186168_arteries_withRoots.tsv", "4332.1132400990/arteries/4332.1132400990_arteries_withRoots.tsv", "4584.1131197859/arteries/4584.1131197859_arteries_withRoots.tsv", "4840.1113050391/arteries/4840.1113050391_arteries_withRoots.tsv")

filenamepath_vec_veins <- c("1396.1132404220/veins/1396.1132404220_veins_withRoots.tsv", "2064.1140250077/veins/2064.1140250077_veins_withRoots.tsv", "2604.1126357612/veins/2604.1126357612_veins_withRoots.tsv", "3300.1120308766/veins/3300.1120308766_veins_withRoots.tsv", "3356.1097932728/veins/3356.1097932728_veins_withRoots.tsv", "3476.1102167727/veins/3476.1102167727_veins_withRoots.tsv", "3840.1105186168/veins/3840.1105186168_veins_withRoots.tsv", "4332.1132400990/veins/4332.1132400990_veins_withRoots.tsv", "4584.1131197859/veins/4584.1131197859_veins_withRoots.tsv", "4840.1113050391/veins/4840.1113050391_veins_withRoots.tsv")

# Inialize lists of dataframes and read into R.
arteries_resorted <- list()
veins_resorted <- list()

for(i in 1:10){
  arteries_resorted[[i]] <- read.csv(file = filenamepath_vec_arteries_resorted[i], header = TRUE)
  veins_resorted[[i]] <- read.csv(file = filenamepath_vec_veins_resorted[i], header = TRUE)
}

# Initialize lists of unique vessel information (removing x,y,z information)
arteries_resorted_unique <- list()
veins_resorted_unique <- list()
resorted_cols <- c("nodeid", "parent", "treeid", "radius_vol", "radius_obs", "length", "volume", "n")

for(i in 1:10){
  temp_arteries <- arteries_resorted[[i]][FALSE,which(names(arteries_resorted[[i]]) %in% resorted_cols)]
  temp_veins <- veins_resorted[[i]][FALSE,which(names(veins_resorted[[i]]) %in% resorted_cols)]
  for(k in unique(arteries_resorted[[i]]$nodeid)){
    temp_arteries <- rbind(temp_arteries, unique(arteries_resorted[[i]][which(arteries_resorted[[i]]$nodeid == k),which(names(arteries_resorted[[i]]) %in% resorted_cols)]))
  }
  for(k in unique(veins_resorted[[i]]$nodeid)){
    temp_veins <- rbind(temp_veins, unique(veins_resorted[[i]][which(veins_resorted[[i]]$nodeid == k),which(names(veins_resorted[[i]]) %in% resorted_cols)]))
  }
  arteries_resorted_unique[[i]] <- temp_arteries
  veins_resorted_unique[[i]] <- temp_veins
  names(arteries_resorted_unique[[i]])[3] <- "radius"
  names(veins_resorted_unique[[i]])[3] <- "radius"
}

rm(arteries_resorted, veins_resorted)

for(i in 1:10){
  names(arteries_resorted_unique[[i]])[3] <- "radius"
  names(veins_resorted_unique[[i]])[3] <- "radius"
}


## Need below block in a loop for each list.  Make new lists, called hs_arteries and hs_veins.
hs_arteries <- list()
hs_veins <- list()

for(i in 1:10){
  # Identify unique tree id values to loop through
  arteries_treeids <- unique(arteries_resorted_unique[[i]]$treeid)
  veins_treeids <- unique(veins_resorted_unique[[i]]$treeid)
  
  # Initialize the final data frame, per list item, by running horton_strahler_relabel on first tree id from treeids vector.
  art_tree <- arteries_resorted_unique[[i]][which(arteries_resorted_unique[[i]]$treeid == arteries_treeids[1]), ]
  vein_tree <- veins_resorted_unique[[i]][which(veins_resorted_unique[[i]]$treeid == veins_treeids[1]), ]

  hs_arteries[[i]] <- horton_strahler_relabel(treeDF = art_tree)
  hs_veins[[i]] <- horton_strahler_relabel(treeDF = vein_tree)
  
  # Loop through remaining treeids.
  for(j in arteries_treeids[-1]){
    art_tree <- arteries_resorted_unique[[i]][which(arteries_resorted_unique[[i]]$treeid == j), ]
    hs_arteries[[i]] <- tryCatch(
      {
        rbind(hs_arteries[[i]], horton_strahler_relabel(treeDF = art_tree))
      },
      error = function(cond) {
        message("nodeid == parent error?")
        message(i)
        message(j)
        message("here is original error message")
        message(cond)
      })
    }
  # Loop through remaining treeids.
  for(j in veins_treeids[-1]){
    vein_tree <- veins_resorted_unique[[i]][which(veins_resorted_unique[[i]]$treeid == j), ]
    hs_veins[[i]] <- tryCatch(
      {
        rbind(hs_veins[[i]], horton_strahler_relabel(treeDF = vein_tree))
      },
      error = function(cond) {
        message("nodeid == parent error?")
        message(i)
        message(j)
        message("here is original error message")
        message(cond)
      })
    }
}

## Now caluclate and add gamma distributions (both regular and horton-strahler)
for(i in 1:length(hs_arteries)){
  gamma <- get.ratio(treeDF = hs_arteries[[i]], treelist = hs_arteries[[i]]$length, is.hs = FALSE)
  hs_gamma <- get.ratio(treeDF = hs_arteries[[i]], treelist = hs_arteries[[i]]$hs_length, is.hs = TRUE)
  hs_arteries[[i]] <- data.frame(hs_arteries[[i]], "gamma" = gamma, "hs_gamma" = hs_gamma)
}

for(i in 1:length(hs_veins)){
  gamma <- get.ratio(treeDF = hs_veins[[i]], treelist = hs_veins[[i]]$length, is.hs = FALSE)
  hs_gamma <- get.ratio(treeDF = hs_veins[[i]], treelist = hs_veins[[i]]$hs_length, is.hs = TRUE)
  hs_veins[[i]] <- data.frame(hs_veins[[i]], "gamma" = gamma, "hs_gamma" = hs_gamma)
}

# build data.frames for all artery and all vein gamma measures.
art_hs_gam_gen <- rbind(hs_arteries[[1]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[2]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[3]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[4]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[5]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[6]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[7]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[8]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[9]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_arteries[[10]][,c("generation", "gen_diff", "gamma", "hs_gamma")])
vein_hs_gam_gen <- rbind(hs_veins[[1]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[2]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[3]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[4]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[5]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[6]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[7]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[8]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[9]][,c("generation", "gen_diff", "gamma", "hs_gamma")], hs_veins[[10]][,c("generation", "gen_diff", "gamma", "hs_gamma")])

#### Reading in and recalculating labelling and scale factors for HHT data. ####

setwd("/Users/alexwork/odrive/Dropbox/Research/plant_animal_vasculature_comparison/datasets/masters/")

hht <- read.csv(file = "hht_master.csv", row.names = 1)
hht$nodeid <- as.character(hht$nodeid)
hht$parent <- as.character(hht$parent)
hht$parent[which(hht$parent == 0)] <- NA
hht <- hht[,c(1:3, 5:ncol(hht))]

hs_hht <- horton_strahler_relabel(treeDF = hht[which(hht$indv == "hht01"),], treeid = "01")

for(i in 2:length(levels(hht$indv))){
  hs_hht <- rbind(hs_hht, horton_strahler_relabel(treeDF = hht[which(hht$indv == levels(hht$indv)[i]),], treeid = strsplit(x = levels(hht$indv)[i], split = "hht")[[1]][[2]]))
}

# hs_hht$hs_nodeid <- as.character(hs_hht$hs_nodeid)
# hs_hht$hs_parent <- as.character(hs_hht$hs_parent)

hs_hht <- data.frame(hs_hht, "hs_gamma" = get.ratio(treeDF = hs_hht, treelist = hs_hht$hs_length, is.hs = TRUE))

#### Reading in and recalculating labelling and scale factors for ML data. ####

ml <- read.csv(file = "mouselung_master.csv", row.names = 1)

ml$nodeid <- as.character(ml$nodeid)
ml$parent <- as.character(ml$parent)
ml$parent[which(ml$parent == 0)] <- NA
ml <- ml[,c(1:2, 4:ncol(ml))]

hs_ml <- horton_strahler_relabel(treeDF = ml, treeid = "01")

hs_ml$hs_nodeid <- as.numeric(hs_ml$hs_nodeid)
hs_ml$hs_parent <- as.numeric(hs_ml$hs_parent)

hs_ml <- data.frame(hs_ml, "hs_gamma" = get.ratio(treeDF = hs_ml, treelist = hs_ml$hs_length, is.hs = TRUE))

#### Reading in and recalculating labelling and scale factors for balsa data. ####

balsa <- read.csv(file = "balsa_master.csv", row.names = 1)

balsa$nodeid <- as.character(balsa$nodeid)
balsa$parent <- as.character(balsa$parent)
balsa$parent[which(balsa$parent == "base")] <- NA
balsa <- balsa[,c(1:2, 4:ncol(balsa))]

hs_balsa <- horton_strahler_relabel(treeDF = balsa, treeid = "01")

hs_balsa$hs_nodeid <- as.numeric(hs_balsa$hs_nodeid)
hs_balsa$hs_parent <- as.numeric(hs_balsa$hs_parent)

hs_balsa <- data.frame(hs_balsa, "hs_gamma" = get.ratio(treeDF = hs_balsa, treelist = hs_balsa$hs_length, is.hs = TRUE))

#### Reading in and recalculating labelling and scale factors for pinon data. ####

pinon <- read.csv(file = "pinon_master.csv", row.names = 1)

pinon$parent[which(pinon$parent == 0.0)] <- NA
pinon <- pinon[,c(1:2, 4:ncol(pinon))]

hs_pinon <- horton_strahler_relabel(treeDF = pinon, treeid = "01")

hs_pinon$hs_nodeid <- as.numeric(hs_pinon$hs_nodeid)
hs_pinon$hs_parent <- as.numeric(hs_pinon$hs_parent)

hs_pinon <- data.frame(hs_pinon, "hs_gamma" = get.ratio(treeDF = hs_pinon, treelist = hs_pinon$hs_length, is.hs = TRUE))


#### Graphing and calculating ####
setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/horton_strahler_relabeling/")

## Lung Arteries
# remove infinities...
art_hs_gam_gen <- art_hs_gam_gen[-c(which(art_hs_gam_gen$gamma == Inf), which(art_hs_gam_gen$hs_gamma == Inf)),]

# plot
png(filename = "human_lung_arteries.png", width = 6, height = 4, units = "in", res = 300)
hist(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], breaks = 800, xlim = c(0, 10), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Human Lung Arteries")
hist(art_hs_gam_gen$gamma, breaks = 1200, xlim = c(0, 10), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(art_hs_gam_gen$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 6, y = 1500, legend = c("HS", paste("median = ", round(x = median(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.696), 3), sep = ""), "Regular", paste("median = ", round(x = median(art_hs_gam_gen$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(1.3), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()

hist(art_hs_gam_gen$gamma, breaks = 500, xlim = c(0, 10))
hist(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], breaks = 400, xlim = c(0, 10))

# means
mean(art_hs_gam_gen$gamma, na.rm = T)
mean(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], na.rm = T)

# medians
median(art_hs_gam_gen$gamma, na.rm = T)
median(art_hs_gam_gen$hs_gamma[which(art_hs_gam_gen$gen_diff == 1)], na.rm = T)

## Lung Veins
# remove infinites
vein_hs_gam_gen <- vein_hs_gam_gen[-c(which(vein_hs_gam_gen$gamma == Inf), which(vein_hs_gam_gen$hs_gamma == Inf)),]
# plot
png(filename = "human_lung_veins.png", width = 6, height = 4, units = "in", res = 300)
hist(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], breaks = 800, xlim = c(0, 10), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Human Lung Veins")
hist(vein_hs_gam_gen$gamma, breaks = 1200, xlim = c(0, 10), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(vein_hs_gam_gen$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 6, y = 1200, legend = c("HS", paste("median = ", round(x = median(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.681), 3), sep = ""), "Regular", paste("median = ", round(x = median(vein_hs_gam_gen$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(1.278), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()


hist(vein_hs_gam_gen$gamma, breaks = 550, xlim = c(0, 10))
hist(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], breaks = 400, xlim = c(0, 10))
#  means
mean(vein_hs_gam_gen$gamma, na.rm = T)
mean(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], na.rm = T)
# medians
median(vein_hs_gam_gen$gamma, na.rm = T)
median(vein_hs_gam_gen$hs_gamma[which(vein_hs_gam_gen$gen_diff == 1)], na.rm = T)

## HHT Arteries
png(filename = "human_head_torso.png", width = 6, height = 4, units = "in", res = 300)
hist(hs_hht$hs_gamma[which(hs_hht$gen_diff == 1)], breaks = 400, xlim = c(0, 5), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Human Head and Torso")
hist(hs_hht$gamma, breaks = 1000, xlim = c(0, 5), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(hs_hht$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(hs_hht$hs_gamma[which(hs_hht$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 3, y = 250, legend = c("HS", paste("median = ", round(x = median(hs_hht$hs_gamma[which(hs_hht$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.435), 3), sep = ""), "Regular", paste("median = ", round(x = median(hs_hht$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.943), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()

hist(hs_hht$gamma, breaks = 900, xlim = c(0, 5))
mean(hs_hht$gamma, na.rm = T)
median(hs_hht$gamma, na.rm = T)

hist(hs_hht$hs_gamma, breaks = 600, xlim = c(0, 5))
mean(hs_hht$hs_gamma, na.rm = T)
median(hs_hht$hs_gamma[which(hs_hht$gen_diff == 1)], na.rm = T)

## ML Arteries
png(filename = "mouse_lung_vessels.png", width = 6, height = 4, units = "in", res = 300)
hist(hs_ml$hs_gamma[which(hs_ml$gen_diff == 1)], breaks = 200, xlim = c(0, 5), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Mouse Lung")
hist(hs_ml$gamma, breaks = 400, xlim = c(0, 5), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(hs_ml$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(hs_ml$hs_gamma[which(hs_ml$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 3, y = 150, legend = c("HS", paste("median = ", round(x = median(hs_ml$hs_gamma[which(hs_ml$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.535), 3), sep = ""), "Regular", paste("median = ", round(x = median(hs_ml$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(1.167), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()

hist(hs_ml$gamma, breaks = 400, xlim = c(0, 10))
mean(hs_ml$gamma, na.rm = T)
median(hs_ml$gamma, na.rm = T)

hist(hs_ml$hs_gamma, breaks = 300, xlim = c(0, 5))
mean(hs_ml$hs_gamma, na.rm = T)
median(hs_ml$hs_gamma[which(hs_ml$gen_diff == 1)], na.rm = T)


## Balsa Branches
png(filename = "balsa_branches.png", width = 6, height = 4, units = "in", res = 300)
hist(hs_balsa$hs_gamma[which(hs_balsa$gen_diff == 1)], breaks = 200, xlim = c(0, 5), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Balsa Tree")
hist(hs_balsa$gamma, breaks = 1000, xlim = c(0, 5), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(hs_balsa$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(hs_balsa$hs_gamma[which(hs_balsa$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 3, y = 140, legend = c("HS", paste("median = ", round(x = median(hs_balsa$hs_gamma[which(hs_balsa$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.301), 3), sep = ""), "Regular", paste("median = ", round(x = median(hs_balsa$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.556), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()


hist(hs_balsa$gamma, breaks = 1000, xlim = c(0, 5))
mean(hs_balsa$gamma, na.rm = T)
median(hs_balsa$gamma, na.rm = T)

hist(hs_balsa$hs_gamma, breaks = 200, xlim = c(0, 5))
mean(hs_balsa$hs_gamma, na.rm = T)
median(hs_balsa$hs_gamma[which(hs_balsa$gen_diff == 1)], na.rm = T)


## Pinon Branches
png(filename = "pinon_branches.png", width = 6, height = 4, units = "in", res = 300)
hist(hs_pinon$hs_gamma[which(hs_pinon$gen_diff == 1)], breaks = 100, xlim = c(0, 5), col = "blue", density = 40, xlab = "Gamma", main = "Histogram of Gamma for Regular and \n Horton-Strahler Labeling, Pinon Tree")
hist(hs_pinon$gamma, breaks = 750, xlim = c(0, 5), add = TRUE, col = "red", density = 20, angle = 135, xlab = "")
abline(v = median(hs_pinon$gamma, na.rm = T), col = "red", lwd = 4)
abline(v = median(hs_pinon$hs_gamma[which(hs_pinon$gen_diff == 1)], na.rm = T), col = "blue", lwd = 4)
legend(x = 3, y = 250, legend = c("HS", paste("median = ", round(x = median(hs_pinon$hs_gamma[which(hs_pinon$gen_diff == 1)], na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.436), 3), sep = ""), "Regular", paste("median = ", round(x = median(hs_pinon$gamma, na.rm = T), digits = 3), sep = ""), paste("D = ", round(-log(2)/log(0.857), 3), sep = "")), fill = c("blue", "blue", "blue",  "red", "red", "red"), density = c(40, 200, 200, 20, 200, 200), angle = c(45, 45, 45, 135, 135, 135))
dev.off()



mean(hs_pinon$gamma, na.rm = T)
median(hs_pinon$gamma, na.rm = T)

mean(hs_pinon$hs_gamma, na.rm = T)
median(hs_pinon$hs_gamma[which(hs_pinon$gen_diff == 1)], na.rm = T)
