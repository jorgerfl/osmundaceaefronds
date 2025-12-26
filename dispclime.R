### *** Quick readme  *** ###

# This script includes functions and code that have been used before
# for computing and plotting a posteriori time-calibrated parsimony trees
# Some functions are utilities for tree comparisons and plotting
# (e.g., *vectoreame* and *gmst_rescaler*), others are for computing the time-calibrated
# trees (*prepare_my_data* and *generate_dated_tree*), and another one for 
# computing and plotting ltt curves along with paleotemperatures (*ltt_and_clime*).
# These functions are kept as legacy/reproducibility and are fully functional.
# 
# The current script adds code to perform disparity-through-time analyses
# (from *Claddis* and *dispRity*). Function *my_htu_dataset* performs ancestral
# character state reconstructions based on an input dated tree and adds HTUs
# to build a matrix with tips and ancestors upon which DTT is calculated (function *process_disp*).

# To replicate analyses, run code lines starting from section 8.

### *** ------------ *** ###

#0000. Load packages ####
# install from CRAN if not installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# install from GitHub if not installed
install_github_if_missing <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github(repo)
  }
}

# CRAN packages
cran_packages <- c(
  "phytools", "ape", "treeio", "TreeTools", "Claddis", "ggtree",
  "tidytree", "tidyverse", "stringr", "ggplot2", "ggnewscale",
  "deeptime", "dplyr", "castor", "paleotree", "svglite", "phangorn",
  "tibble", "ggplotify", "tidyr", "plyr", "scales", "patchwork",
  "zoo"
)

# Install CRAN packages if missing
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# GitHub-only packages
# (dispRity is not on CRAN, must be installed from GitHub)
install_github_if_missing("dispRity", "TGuillerme/dispRity")

# Finally, load everything
libs <- c(cran_packages, "dispRity")
invisible(lapply(libs, library, character.only = TRUE))


# library(phytools)
# library(ape)
# library(treeio)
# library(TreeTools)
# library(Claddis)
# library(ggtree)
# library(tidytree)
# library(tidyverse)
# library(stringr)
# library(ggplot2)
# library(ggnewscale)
# library(deeptime)
# library(dplyr)
# library(castor)
# library(paleotree)
# library(svglite)
# library(phangorn)
# library(tibble)
# library(ggplotify)
# library(tidyr)
# library(plyr)
# library(scales)
# library(patchwork)
# library(zoo)
# library(dispRity)



#000. FUNCTIONS ####

# Tree-comparison and extraction utility.
# After matching nodes and tips, extracts node ages and branch lengths and returns a matrix
vectoreame <- function(dated_multiphylo_list, dated_multiphylo_tibble_list, reference_phylo, reference_tibble, A){
  
  if(!inherits(dated_multiphylo_list, "multiPhylo")){ stop("Not multiPhylo list object given")}
  if(!is.list(dated_multiphylo_tibble_list)){ stop("Not dated trees tibble list object given")}
  
  if(!inherits(reference_phylo, "phylo")){ stop("Not reference Phylo object given")}
  if(!inherits(reference_tibble, "tbl_tree")){ stop("Not reference tibble tree object given")}
  
  #Actual code starts here
  eqterminals <- vector(mode = "logical", length = length(reference_phylo$tip.label))
  eqterminals[] <- NA
  tr2labels <- matchLabels(reference_phylo, dated_multiphylo_list[[A]])
  tr2labels <- as.vector(tr2labels[,2])
  eqterminals[1:length(tr2labels)] <- tr2labels
  
  eqnodes <- vector(mode = "logical", length = (nrow(reference_tibble)-length(reference_phylo$tip.label)) )
  eqnodes[] <- NA
  tr2nodes <- matchNodes(reference_phylo, dated_multiphylo_list[[A]], method = "descendants")
  tr2nodes <- as.vector(tr2nodes[,2])
  eqnodes[1:length(tr2nodes)] <- tr2nodes
  
  eqgroups <- c(eqterminals,eqnodes)
  eqages <- vector(mode = "integer", length = length(eqgroups))
  eqages[] <- NA
  
  eqbranchlen <- vector(mode = "integer", length = length(eqgroups))
  eqbranchlen[] <- NA
  
  
  for (zz in 1:length(eqgroups)) {
    if(!is.na(eqgroups[zz])){
      eqages[zz] <- dated_multiphylo_tibble_list[[A]]$nodeages[eqgroups[zz]]
      eqbranchlen[zz] <- dated_multiphylo_tibble_list[[A]]$branch.length[eqgroups[zz]]
    }
  }
  #end of actual code
  
  stringA <- paste("node", A, sep = ".")
  stringB <- paste("ages", A, sep = ".")
  stringC <- paste("branch.length", A, sep = ".")
  
  output_matrix <- matrix(data = c(eqgroups, eqages, eqbranchlen), ncol = 3)
  colnames(output_matrix) <- c(stringA, stringB, stringC)
  
  return(output_matrix)
}

#To process LTT and paleoclimate data. It uses Judd et al climate data
#As compared to traditional LTT plots, this function "smooths" the LTT profile
ltt_and_clime <- function(mytrees, clime="PhanDA_GMSTandCO2_percentiles.csv"){
  
    #Output: 
  #A list with:
  #'lttclime' tibble that stores median ltt and median temperature, and their respective time points.
  #'ltt_ribbon_data' dataframe with 95% limits and time points for creating a ribbon around median ltt
  #'clim_ribbon_data' dataframe with 95% limits and time points for creating a ribbon around median temperature   
  
  flag_multiphylo = 0
  flag_phylo = 0
  
  #LTT from mytrees
  if(inherits(mytrees, "multiPhylo")){
    cat("Handling more than a single input tree - they are assumed to be the same tree but with different inferred ages\n")
    flag_multiphylo = 1
    mytree_ltt <- ltt.plot.coords(mytrees[[1]])
    mytree_ltt[,1] <- mytree_ltt[,1] * (-1)
    colnames(mytree_ltt)[1] <- paste("_","1",".", colnames(mytree_ltt)[1], sep = "")
    colnames(mytree_ltt)[2] <- paste("_","1",".", colnames(mytree_ltt)[2], sep = "")
    for (i in 2:length(mytrees)) {
      toadd <- ltt.plot.coords(mytrees[[i]])
      toadd[,1] <- toadd[,1] * (-1)
      colnames(toadd)[1] <- paste("_",i,".", colnames(toadd)[1], sep = "")
      colnames(toadd)[2] <- paste("_",i,".", colnames(toadd)[2], sep = "")
      mytree_ltt <- cbind(mytree_ltt, toadd)
    }
    
    mytree_ltt <- as.data.frame(mytree_ltt)
    
    #Select LTT columns
    median_ltt <- mytree_ltt %>% select(contains(".N"))
    #Calculate median + 95% CI (2.5% and 97.5% quantiles) for each row
    ltt.row_stats <- t(apply(median_ltt, 1, function(x) {
      quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
    }))
    #Rename columns for clarity
    colnames(ltt.row_stats) <- c("median.ltt", "ltt.lower_95", "ltt.upper_95")
    
    #Select time columns
    median_time <- mytree_ltt %>% select(contains(".time"))
    #Calculate median + 95% CI (2.5% and 97.5% quantiles) for each row
    time.row_stats <- t(apply(median_time, 1, function(x) {
      quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
    }))
    #Rename columns for clarity
    colnames(time.row_stats) <- c("median.time", "time.lower_95", "time.upper_95")
    
    summary_ltt <- cbind(ltt.row_stats, time.row_stats)
    summary_ltt <- as.data.frame(summary_ltt)
    
    # Create a separate dataset for geom_ribbon
    ltt_ribbon_data <- summary_ltt %>%
      select(median.time, ltt.lower_95, median.ltt, ltt.upper_95) %>%
      drop_na()
    # Where the lower limits have no difference with median, I create an artificiial lower limit of 1
    difflow = ltt_ribbon_data[[3]] - ltt_ribbon_data[[2]]
    difflow_index = which(difflow == 0)
    ltt_ribbon_data[[2]][difflow_index] <- ltt_ribbon_data[[3]][difflow_index] - 1
    # Where the upper limits have no difference with median, I create an artificiial lower limit of 1
    diffup = ltt_ribbon_data[[4]] - ltt_ribbon_data[[3]]
    diffup_index = which(diffup == 0)
    ltt_ribbon_data[[4]][diffup_index] <- ltt_ribbon_data[[3]][diffup_index] + 1
    
    
  } #end of multiPhylo
  else{
    cat("There is a single input tree\n")
    flag_phylo = 1
    mytree_ltt <- ltt.plot.coords(mytrees)
    mytree_ltt[,1] <- mytree_ltt[,1] * (-1)
    colnames(mytree_ltt)[1] <- paste("median",".", colnames(mytree_ltt)[1], sep = "")
    colnames(mytree_ltt)[2] <- paste("median.", colnames(mytree_ltt)[2], sep = "")
    
    colnames(mytree_ltt)[1] <- "median.time"
    colnames(mytree_ltt)[2] <- "median.ltt"
    
    summary_ltt <- mytree_ltt
  }
  
  #Paleoclimate dataset
  
  # Read paleoclimate data from Judd et al. (2024): 10.1126/science.adk3705
  paleodata <- read.csv(clime, header = TRUE, sep = ",")
  paleodata <- paleodata[, (names(paleodata) %in% c("AverageAge", "GMST_05", "GMST_50", "GMST_95"))]
  colnames(paleodata)[1] <- "clim_time"
  
  cat(colnames(paleodata))
  
  # Create a separate dataset for geom_ribbon
  paleo_ribbon_data <- paleodata %>%
    select(clim_time, GMST_05, GMST_50, GMST_95) %>%
    drop_na()
  
  # Combining with summary_ltt
  combined_paleodata_ltt <- merge(summary_ltt, paleodata, all = TRUE) %>% 
    drop_na()
  
  # Reshape data for being plotted
  paleoltt_long_data <- combined_paleodata_ltt %>%
    pivot_longer(
      cols = c(median.ltt, GMST_50),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(time = case_when(
      variable == "median.ltt" ~ median.time,
      variable == "GMST_50" ~ clim_time
    ))
  
  
  if(flag_multiphylo == 1) {
    return(list(lttclime = paleoltt_long_data, clim_ribbon_data = paleo_ribbon_data, ltt_ribbon_data = ltt_ribbon_data))
  } else if (flag_phylo == 1)
  {return(list(lttclime = paleoltt_long_data, clim_ribbon_data = paleo_ribbon_data, ltt_ribbon_data = NULL))}

}

#Function to rescale GMST values relative to LTT values
gmst_rescaler <- function(ltt_data, gmst_data) {
  function(x) {
    (x - min(ltt_data$value, na.rm = TRUE)) /
      (max(ltt_data$value, na.rm = TRUE) - min(ltt_data$value, na.rm = TRUE)) *
      (max(gmst_data$GMST_50, na.rm = TRUE) - min(gmst_data$GMST_50, na.rm = TRUE)) +
      min(gmst_data$GMST_50, na.rm = TRUE)
  }
}

#Function to prepare data
prepare_my_data <- function(optimal_tree_file, suboptimal_tree_file, FAD_file = "FAD.csv", age_filter = 6) {
  
  # Load optimal and suboptimal trees
  optimal_tree <- ReadTntTree(optimal_tree_file)
  suboptimal_trees <- ReadTntTree(suboptimal_tree_file)
  
  # Sample suboptimal trees if there are 10 or more
  if (inherits(suboptimal_trees, "multiPhylo") && length(suboptimal_trees) >= 10) {
    suboptimal_trees <- sample(suboptimal_trees, 10, replace = FALSE)
  }
  
  # checks whether I have one or multiple trees in "optimal_tree"
  if(inherits(optimal_tree, "multiPhylo")) {
    cat("The object is of class 'multiPhylo' -- a strict consensus is computed \n")
    
    #Establish a consensus tree, to function as reference
    reference_tree <- consensus(optimal_tree, rooted = TRUE)
    
  }
  else if(inherits(optimal_tree, "phylo")){
    cat("The object is of class 'phylo' -- no strict consensus computed \n")
    reference_tree <- optimal_tree
  }
  else
  { stop("The object is neither 'phylo' nor 'multiPhylo'.\n") }
  
  # all individual trees
  allreplicates <- c(optimal_tree, suboptimal_trees)
  
  
  
  # Process the FAD-LAD matrix
  FAD <- read.csv(FAD_file, header = TRUE, row.names = 1)
  matched_indices <- match(reference_tree$tip.label, rownames(FAD))
  valid_indices <- !is.na(matched_indices)
  FAD_filtered <- FAD[matched_indices[valid_indices], ]
  
  # Apply the FAD/LAD filter using the new argument
  FAD_filtered$FAD[FAD_filtered$FAD < age_filter] <- 0
  FAD_filtered$LAD[FAD_filtered$LAD < age_filter] <- 0
  
  # Convert to matrix
  rangeCont <- as.matrix(FAD_filtered)
  
  # Return all three outputs as a list
  return(list(allreplicates = allreplicates, reference_tree = reference_tree, optimal_tree = optimal_tree, rangeCont = rangeCont))
}

#Function to date trees (it calls other functions)
generate_dated_tree <- function(processed_data = my_processed_data, tree_tibble = reference_tree_tibble, bmass = TRUE, others = NA){
  
  dated_trees <- NULL
  for (i in 1:length(processed_data$allreplicates)) {
    tree_i <- processed_data$allreplicates[[i]]
    tree_i_dated <- timePaleoPhy(tree = tree_i,
                                 ntrees = sampletime, 
                                 randres = F, 
                                 vartime = vartimer,
                                 dateTreatment = "minMax",
                                 timeData = processed_data$rangeCont, 
                                 type = "equal", 
                                 plot = F)
    
    dated_trees <- append(dated_trees, tree_i_dated, after = length(dated_trees)) }
  
  #Create a list of tibbles, where each tibble is a dated tree
  dated_trees_tibbled <- vector(mode = "list", length = length(dated_trees))
  for (z in 1:length(dated_trees_tibbled)) {
    dated_trees_tibbled[[z]] <- as_tibble(dated_trees[[z]])
    dated_trees_tibbled[[z]][["nodeages"]] <- dateNodes(tree = dated_trees[[z]])
    dated_trees_tibbled[[z]][["nodeages"]] <- round(dated_trees_tibbled[[z]][["nodeages"]], digits = 3)
    dated_trees_tibbled[[z]][["nodeages"]] <- as.vector(dated_trees_tibbled[[z]][["nodeages"]]) }
  
  #Now time-scale only the optimal trees to compute consensus edge.lengths in the reference tree
  optimal_dated <- NULL
  if(inherits(processed_data$optimal_tree, "phylo")){
    
    tree_i_optimal <- processed_data$optimal_tree
    tree_i_dated_optimal <- timePaleoPhy(tree = tree_i_optimal,
                                         ntrees = sampletime, 
                                         randres = F, 
                                         vartime = vartimer,
                                         dateTreatment = "minMax",
                                         timeData = processed_data$rangeCont, 
                                         type = "equal", 
                                         plot = F)
    optimal_dated <- append(optimal_dated, tree_i_dated_optimal, after = length(optimal_dated))
    
  } else if(inherits(processed_data$optimal_tree, "multiPhylo")) {
    
    for (i in 1:length(processed_data$optimal_tree)) {
      tree_i_optimal <- processed_data$optimal_tree[[i]]
      tree_i_dated_optimal <- timePaleoPhy(tree = tree_i_optimal,
                                           ntrees = sampletime, 
                                           randres = F, 
                                           vartime = vartimer,
                                           dateTreatment = "minMax",
                                           timeData = processed_data$rangeCont, 
                                           type = "equal", 
                                           plot = F)
      optimal_dated <- append(optimal_dated, tree_i_dated_optimal, after = length(optimal_dated))} 
  }
  consedges <- consensus.edges(optimal_dated)
  
  
  final_mat <- NULL
  for (j in 1:length(dated_trees)) {
    temporal_mat <- vectoreame(dated_multiphylo_list = dated_trees, dated_multiphylo_tibble = dated_trees_tibbled,
                               reference_phylo = processed_data$reference_tree, reference_tibble = tree_tibble, A = j)
    final_mat <- cbind(final_mat, temporal_mat) }
  
  # Remove columns with names matching "node.X" and "ages.X"
  final_branch.lengths <- final_mat[, !grepl("^node\\..*$", colnames(final_mat))]
  final_branch.lengths <- final_branch.lengths[, !grepl("^ages\\..*$", colnames(final_branch.lengths))]
  
  # Estimate median of branch lengths as final columns
  final_branch.lengths <- cbind(final_branch.lengths, NA)
  colnames(final_branch.lengths)[ncol(final_branch.lengths)] <- "branch.length"
  for (yy in 1:nrow(final_branch.lengths)) {
    final_branch.lengths[yy,ncol(final_branch.lengths)] <- median(unname(unlist(final_branch.lengths[yy,1:ncol(final_branch.lengths)])),na.rm = T) 
  }
  
  tree_tibble <- cbind(tree_tibble, final_mat)
  
  # Remove columns with names matching "node.X" and "branch.length.X"
  tree_tibble <- tree_tibble[, !grepl("^node\\..*$", colnames(tree_tibble))]
  tree_tibble <- tree_tibble[, !grepl("^branch\\.length\\..*$", colnames(tree_tibble))]
  
  # Estimate median values from ages, minimum and maximum ages and median branch lengths
  for (ii in 1:nrow(tree_tibble)) {
    
    #Ages:
    tree_tibble$median_age[ii] <- median(unname(unlist(tree_tibble[ii,8:ncol(tree_tibble)])),
                                         na.rm = T)
    tree_tibble$min_age[ii] <- min(unname(unlist(tree_tibble[ii,8:ncol(tree_tibble)])),
                                   na.rm = T)
    tree_tibble$max_age[ii] <- max(unname(unlist(tree_tibble[ii,8:ncol(tree_tibble)])),
                                   na.rm = T)  }
  
  
  # Final plotting preparation
  main_tibble <- tree_tibble[,1:6]
  main_tibble <- cbind(main_tibble,final_branch.lengths[,ncol(final_branch.lengths)])
  colnames(main_tibble)[ncol(main_tibble)] <- "branch.length"
  
  #Create ranges from min and max ages and exclude them
  range_string <-  list()
  for (x in 1:nrow(main_tibble)) {
    
    range_string[[x]] <- c(main_tibble$min_age[[x]], main_tibble$max_age[[x]])
    
  }
  
  cols_to_delete <- c("min_age", "max_age")
  main_tibble <- main_tibble[,!colnames(main_tibble) %in% cols_to_delete]
  main_tibble$range <- range_string
  
  main_tibble <- as_tibble(main_tibble)
  write.nexus(consedges, file = "temporal.temp")
  main_tree <- as_tibble(TNTOrder(read.nexus("temporal.temp")))
  main_tree$median_ages <- main_tibble$median_age
  main_tree$range <- main_tibble$range
  for (i in 1:length(main_tree$label)) { 
    if(main_tree$label[i] == "1") {
      main_tree$label[i] <- paste("HTU_", i, sep = "")
    }
  }
  
  
  #This lines are for including body mass reconstructions from TNT
  #bmass is taken from a tab delimited file with taxa, min body mass, max body mass, and htu body mass values
  #if you wish to include other data, follow the same structure
  if(isTRUE(bmass)){
    
    bm <- read.table("bm.csv", quote="\"", comment.char="")
    bm$V2 <- as.numeric(bm$V2)
    bm$V3 <- as.numeric(bm$V3)
    for (j in 1:nrow(bm)) { if(grepl("HTU_", bm$V1[j])){ bm$V1[j] <- paste("HTU_", j, sep = "") } }
    for (j in 1:nrow(bm)) { bm$median_body_mass[j] <- median(c(bm$V2[j],bm$V3[j])) }
    
    main_treerepl <- as.treedata(main_tree)
    bmtnt <- as_tibble(TNTOrder(main_treerepl@phylo))
    bm_reorder <- bm[match(bmtnt$label, bm$V1), ] #body mass dataset from TNT is reorder based on the tibble
    main_tree$body_mass <- as.numeric(bm_reorder$median_body_mass)
    #for (j in 51:70) { main_tree$body_mass[j] <- NA }
  } #End of lines for including body mass reconstructions from TNT
  
  
  #As with boddy mass, but adding other continuous trait as stored in the given file
  if(!is.na(others)) {
    if(file.exists(others)) {
      
      other_data <- read.table(others, quote="\"", comment.char="")
      other_data$V2 <- as.numeric(other_data$V2)
      other_data$V3 <- as.numeric(other_data$V3)
      for (j in 1:nrow(other_data)) { if(grepl("HTU_", other_data$V1[j])){ other_data$V1[j] <- paste("HTU_", j, sep = "") } }
      for (j in 1:nrow(other_data)) { other_data$median_trait[j] <- median(c(other_data$V2[j],other_data$V3[j])) }
      
      main_treerepl <- as.treedata(main_tree)
      othertnt <- as_tibble(TNTOrder(main_treerepl@phylo))
      other_data_reorder <- other_data[match(othertnt$label, other_data$V1), ] #dataset is reorder based on the tibble
      main_tree$median_trait <- as.numeric(other_data_reorder$median_trait)
      #for (j in 51:70) { main_tree$median_trait[j] <- NA }
      
      # Process the data
    } else {
      warning(paste("File not found:", others))
    }
  } #End of the lines for including another continuous trait
  
  main_tree <- as.treedata(main_tree)
  
  if (isTRUE(bmass)) {
    return(list(main_tree = main_tree, dated_trees = dated_trees, bmtnt = bmtnt))
  } else {
    return(list(main_tree = main_tree, dated_trees = dated_trees, bmtnt = NA))
  }
  
}

#Creates a dataset with htus and tips, using a dated tree
my_htu_dataset <- function(input_tree, input_data, from_file=FALSE){
  
  
  # read matrix (no htus) from a read_nexus_matrix output
  my_matrix <- input_data
  
  # Tree is to be read from file? Well, I will also need a root age and a file name
  if(isTRUE(from_file) ){
    R_tree_all <- read.mrbayes(file = input_tree)
    R_tree <- R_tree_all@phylo
    R_tree$root.time <- max(as.integer(R_tree_all@data$age_median)) 
  } else {
    R_tree_all <- input_tree #This is a treedata object and should have branch lengths (I am not checking it)
    R_tree <- R_tree_all@phylo
    R_tree$root.time <- max(as.integer(R_tree_all@data$age_median)) 
  }
  
  
  #---Ancestral reconstruction per character and ordination---#
  htu_data <- estimate_ancestral_states(cladistic_matrix = my_matrix, time_tree = R_tree)
  htu_data_ordinated <- Claddis.ordination(htu_data, add = T)
  R_tree$node.label <- row.names(htu_data_ordinated)[(length(R_tree$tip.label)+1):(nrow(htu_data_ordinated))]
  #---End of reconstruction---# 
  
  all <- list("htu_data" = htu_data, "htu_ordinated_data" = htu_data_ordinated, "R_tree" = R_tree)
  
  #return dataset with added htu, and the tree
  return(all)
}

#Estimate disparity using dispRity upon a ordinated matrix (using Claddis)
process_disp <- function(input, time_split = NULL){
  
  #time slicing
  time_slice <- chrono.subsets(
    data = input$htu_ordinated_data, 
    tree = input$R_tree,
    inc.nodes = T,
    t0 = F,
    method = "continuous",
    model = "proximity",
    time = seq(from = input$R_tree$root.time, to = 0, by = (-20))
    
  )
  
  #bootsrapping
  time_sliceBS <- boot.matrix(
    data = time_slice, 
    bootstraps = 1000,
    boot.type = "full")
  
  #estimate dtt
  disparA <- dispRity(data = time_sliceBS, 
                      metric = c(mean, span.tree.length)) #density

  # report <- summary(disparA)
  # mylims <- vector(mode = "integer", length = length(myages))
  # for (i in seq_along(myages)) {
  #   ratio_diff <- abs(( as.numeric(report$subsets) / myages[i]) - 1)
  #   mylims[i] <- which.min(ratio_diff)
  # }
  
  
  disparB <- dispRity(data = time_sliceBS, 
                      metric = c(mean, displacements),
                      reference=c(0, 1)) #position
  
  # for(a in 1:length(disparB$disparity)) {
  #   
  #   position_val <- disparB$disparity[[a]][[2]][1:ncol(disparB$disparity[[a]][[2]])]
  #   position_val <- log(position_val)
  #   disparB$disparity[[a]][[2]] <- position_val
  #   
  # }
  
  
  disparC <- dispRity(data = time_sliceBS, 
                      metric = c(sum, ranges)) #size
  
  metrics <- list("density" = disparA, "position" = disparB, "size" = disparC)
  return(metrics)
}

#Partition cladistic matrix according to type (binary vs multistate) or custom
splitme <- function(mydata, bytype=T, chars = NULL){
  
  if(isTRUE(bytype)){
    binarychars <- vector(mode = "numeric")
    allchars <- seq(1:ncol(mydata[["matrix_1"]][["matrix"]]))
    
    for (ch in 1:ncol(mydata[["matrix_1"]][["maximum_values"]])) {
      chstate_max <- max(as.numeric(mydata[["matrix_1"]][["maximum_values"]][ch]), na.rm = T)
      if(chstate_max <= 1) {
        binarychars <- append(binarychars, ch, after = length(binarychars))
      }
    }
   non_binarychars_index <-  match(binarychars, allchars)
   non_binarychars <- allchars[-non_binarychars_index]
   
   #Now generate three matrices: entire, paetition1-only, and partition2-only
   entire <- mydata
   partition_bin <- mydata
   partition_mult <- mydata
   
   partition_bin[["matrix_1"]][["matrix"]]  <- mydata[["matrix_1"]][["matrix"]][,binarychars]
   partition_bin[["matrix_1"]][["ordering"]]  <- mydata[["matrix_1"]][["ordering"]][binarychars] 
   partition_mult[["matrix_1"]][["matrix"]] <- mydata[["matrix_1"]][["matrix"]][,non_binarychars]
   partition_mult[["matrix_1"]][["ordering"]]  <- mydata[["matrix_1"]][["ordering"]][non_binarychars] 

   return(list("entire" = entire, "bin" = partition_bin,"nonbin"= partition_mult))
   
  } else if (!is.null(chars)){
    entire <- mydata
    partition_custom <- mydata
    partition_custom[["matrix_1"]][["matrix"]]  <- mydata[["matrix_1"]][["matrix"]][,chars]
    partition_custom[["matrix_1"]][["ordering"]]  <- mydata[["matrix_1"]][["ordering"]][chars] 
    partition_custom[["matrix_1"]][["character_weights"]]  <- mydata[["matrix_1"]][["character_weights"]][chars]
    partition_custom[["matrix_1"]][["minimum_values"]]  <- mydata[["matrix_1"]][["minimum_values"]][chars]
    partition_custom[["matrix_1"]][["maximum_values"]]  <- mydata[["matrix_1"]][["maximum_values"]][chars]
    
    return(list("entire" = entire,"custom"= partition_custom))
  }

}

#Plot using plot.dispRity
plot_my_plate <- function(mydisp_object, title=NULL, hformat=F, sqformat=T, vformat=F){
  
  # Plots from the entire matrix
  # Save original graphical settings
  old_par <- par(no.readonly = TRUE)
  
  if(isTRUE(sqformat)){
    par(
      mfrow = c(2, 2),    # Use mfrow for simplicity
      mar = c(5, 5, 2, 1), # Inner margins: bottom, left, top, right
      oma = c(1.5, 1.75, 1.75, 1.5), # Outer margins: gives the whole grid some space
      cex.axis = 1,      # Shrink axis number font size
      cex.lab = 1.75,       # Shrink axis label font size
      bty = "n"            # No box around plots
    )
    
  } else if(isTRUE(hformat)) {
    par(
      mfrow = c(1, 3),    # Use mfrow for simplicity
      mar = c(5, 4, 2, 1), # Inner margins: bottom, left, top, right
      oma = c(1.5, 1.75, 1.75, 1.5), # Outer margins: gives the whole grid some space
      cex.axis = 1,      # Shrink axis number font size
      cex.lab = 1.75,         # Shrink axis label font size
      bty = "n"            # No box around plots
    )
  } else if(isTRUE(vformat)) {
    par(
      mfrow = c(3, 1),    # Use mfrow for simplicity
      mar = c(5, 4, 2, 1), # Inner margins: bottom, left, top, right
      oma = c(1.5, 1.75, 1.75, 1.5), # Outer margins: gives the whole grid some space
      cex.axis = 1,      # Shrink axis number font size
      cex.lab = 1.75,         # Shrink axis label font size
      bty = "n"            # No box around plots
    )
    
  }
  
  # Your plots
  plot.dispRity(mydisp_object$density, type = "continuous", col = c("red", "lightpink", "pink"), elements = F) ; title(main = "Density", col.main = "black", cex.main = 1.75, font.main = 1)
  plot.dispRity(mydisp_object$position, type = "continuous", col = c("red", "lightpink", "pink"), elements = F) ; title(main = "Position", col.main = "black", cex.main = 1.75, font.main = 1)
  plot.dispRity(mydisp_object$size, type = "continuous", col = c("red", "lightpink", "pink"), elements = F) ; title(main = "Size", col.main = "black", cex.main = 1.75, font.main = 1)
  
  # You can add a main title in the outer margin
  mtext(title, outer = TRUE, cex = 1.0, line = 0.5)
  
  # Restore the original settings
  par(old_par)
  
}

#EXTRACT °C FROM PALEODATA
extract_temp <- function(clime = "PhanDA_GMSTandCO2_percentiles.csv", limits = NA, split = NA, root.time = NA){
  
  # Read paleoclimate data from Judd et al. (2024): 10.1126/science.adk3705
  paleodata <- read.csv(clime, header = TRUE, sep = ",")
  paleodata <- paleodata[, (names(paleodata) %in% c("AverageAge", "GMST_05", "GMST_50", "GMST_95"))]
  colnames(paleodata)[1] <- "clim_time"
  paleodata$clim_time <- round(paleodata$clim_time, digits = 2)
  
  cat("The minimum and maximum time ages available are:", min(paleodata$clim_time), max(paleodata$clim_time), "\n")
  
  low <- min(limits)
  upper <- max(limits)
  
  # Subset data to include only rows within the specified time limits
  # limits[1] = lower time bound, limits[2] = upper time bound
  paleodata_subset <- paleodata[paleodata$clim_time >= low & 
                                  paleodata$clim_time <= upper, ]
  
  cat("The ages you picked up range from ", min(paleodata_subset$clim_time), "to ", max(paleodata_subset$clim_time), "\n")
  
  cat("var <- v(")
  for (n in 1:(nrow(paleodata_subset)-1)) {
    cat(paleodata_subset$GMST_50[n])
    cat(", ")
  }
  cat(paleodata_subset$GMST_50[nrow(paleodata_subset)],")\n")
  cat("MAX_VAR_AGE = ", paleodata_subset$clim_time[nrow(paleodata_subset)], "\n")
  
  all <- list("paleo" = paleodata_subset, "limits" = limits)
  return(all)
}

#PLOT SEVERAL TREES
plot_multiple_mrbayes_trees <- function(myinput_trees, ncol = 2, nrow = NULL, titles = NULL, filename = NULL, nexfile=F) {
  # Create empty list
  empty_list <- list()
  
  # Loop through all trees (FIX: 1:length() instead of just length())
  for(i in 1:length(myinput_trees)){
    # Read the MrBayes tree
    if(isTRUE(nexfile))
       { tree <- read.nexus(myinput_trees[i]) }
     else
       { tree <- read.mrbayes(myinput_trees[i]) }
    
    # Create ggtree plot
    treeplot <- ggtree(tree, 
                       layout = "rectangular", 
                       ladderize = TRUE, 
                       right = TRUE,
                       branch.length = "none", 
                       size = 0.05) + 
      geom_tiplab(size = 2, 
                  color = "black",  
                  offset = 0.1) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "mm")) +  # Reduce margins
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.38)))
    # Add to list
    empty_list[[i]] <- treeplot
  }
  
  # Set default titles if not provided
  if (is.null(titles)) {
    titles <- basename(myinput_trees)  # Use filenames as titles
  }
  
  # Add titles to plots
  for (i in 1:length(empty_list)) {
    empty_list[[i]] <- empty_list[[i]] + 
      ggtitle(titles[i]) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  }
  
  # Calculate layout if nrow not specified
  if (is.null(nrow)) {
    nrow <- ceiling(length(empty_list) / ncol)
  }
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(empty_list, 
                              ncol = ncol, 
                              nrow = nrow)
  
  # Print to screen
  print(combined_plot)
  
  # Save to file if filename provided
  if (!is.null(filename)) {
    # Calculate appropriate dimensions
    plot_width <- ncol * 4
    plot_height <- nrow * 3
    ggsave(filename = filename,
           plot = combined_plot,
           width = plot_width,
           height = plot_height,
           dpi = 300)
    message(paste("Plot saved as:", filename))
  }
  
  return(invisible(combined_plot))
}


###* ************************** BLOCK 1 ************************** ####
#*   ******************** A POSTERIORI TIME CALIBRATE **********     #
###* ************************************************************  ###


#0.  STEP - INITIALIZING SOME VARIABLES FOR TIME SCALING ####
input_A = "5dep.tre"
input_B = "suboptimalDEP_5.tre"
vartimer = 10
sampletime = 1000

# 1. Create the data frame for the rectangles
rect_data <- data.frame(
  xmin = c(-193, -49),
  xmax = c(-175, -16),
  ymin = -Inf,
  ymax = Inf
)
my_alpha = 0.55
inferior_limit = (-312)
breaks_vec <- c(inferior_limit,
                -193, -175, #Pleinsbachian-Toarcian, Toarcian-Aalenian (PTo-E and T-OAE events)
                -49, -16    #Opening of the Drake Passage and thermal isolation of Antarctica
                #Even though its opening it is assumed to have concluded by 16-17 Ma,
                #its beginning is blurry but estimated to be around the Eocene-Oligocene transition.
                #See:
                #Caruthers et al. (2013):https://doi.org/10.1016/j.palaeo.2013.05.010
                #Scher & Martin (2006): 10.1126/science.1120044
                #Vincze et al. (2021): https://doi.org/10.1038/s41598-021-99123-0
)


#rgbfill <- rgb(red = 100, green = 200, blue = 250, maxColorValue = 255)
rgbfill <- rgb(red = 120, green = 180, blue = 220, maxColorValue = 255)

# 1. Create the data frame for the rectangles

#1.  STEP - PREPARING DATA ####
my_processed_data <- prepare_my_data(input_A, input_B)


#Create a tibble from reference and plot
reference_tree_tibble <- as_tibble(my_processed_data$reference_tree)
reference_tree_tibble$median_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
reference_tree_tibble$min_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
reference_tree_tibble$max_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
reference_tree_tibble$branch.length <- vector(mode = "logical", length = nrow(reference_tree_tibble))
ggtree(as.phylo(reference_tree_tibble), branch.length = "none") + geom_tiplab(size=2)


#2.  STEP - TIME-SCALING INDIVIDUAL TREES ####
my_dated_data <- generate_dated_tree(others = "rga.csv")
main_tree <- my_dated_data$main_tree

#2B. STEP - PLOTTING TREES ####

main_tree@data$age_mean <- as.double(main_tree@data$age_mean)
main_tree@data$age_median <- as.double(main_tree@data$age_median)

# Replace "_" by blank
main_tree@phylo$tip.label <- gsub("_", " ", main_tree@phylo$tip.label)

#Preserve data for nodes of interest only (those included in node 70)
ingroup_index = unlist(Descendants(main_tree@phylo, 70, type = "all"))
main_tree@data$median_trait[-ingroup_index] <- NA
main_tree@data$body_mass[-ingroup_index] <- NA
main_tree@data$`rateTK02Brlens{1}_median` <- as.double(main_tree@data$`rateTK02Brlens{1}_median`)

# ggtree objects have a data slot that is a data frame.
tree_data_filtered <- ggtree(main_tree)$data |> dplyr::filter(!isTip)
tree_data_filtered$age_median <- as.integer(tree_data_filtered$age_median)

# Start the plot with ggtree
f3 <- ggtree(main_tree, ladderize = TRUE, right = TRUE) +
  theme_tree2() +
  # Add the geom_rect layer first
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +
  # Add the tree and other tree-related layers next
  geom_tree(layout = "rectangular", size = 0.95, aes(colour=`rateTK02Brlens{1}_median`))+ #colour = rgb(99, 99, 99, maxColorValue = 255)) +
  #labs(colour = "Body mass") +
  geom_text(aes(x = branch, label = ifelse(round(age_median, 2) > 0, round(age_median, 2), "")), size = 2.75, hjust=1, vjust=1.3)+
  # Add all other geoms
  # Pass the pre-filtered data frame directly to the data argument
  geom_range(data = tree_data_filtered, range = 'age_0.95HPD', colour = 'red', size = 1.25, alpha = 0.3,
             position = position_nudge(x = -0.5)) +
  scale_color_continuous(low = 'blue', high = 'darkred') +
  #scale_color_viridis_b(option = "magma")+
  geom_tiplab(size = 3.0, offset = 0.5, fontface = "italic") +
  # Finally, add the coordinate and scale transformations
  coord_geo(xlim = c(inferior_limit, Ntip(main_tree) + 1), ylim = c(0.5, Ntip(main_tree)), expand = TRUE,
            dat = list("eras", "periods", "epochs"), abbrv = list(F, F, T),
            skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
            pos = list("top", "top", "top"), alpha = 1, height = unit(0.75, "line"),
            size = "auto", fittext_args = list(size = 22), rot = 0, neg = TRUE) +
  scale_x_continuous(breaks = breaks_vec, labels = abs(breaks_vec),
                     expand = expansion(mult = c(0.009, 0.14))); f3_final <- revts(f3); f3_final

ggsave(filename = "plotted_tree.svg", device = "svg", width = 40, height = 25, units = "cm")


aposttree <- main_tree
main_tree@phylo$node.label <- NULL
write.tree(main_tree@phylo, file = "exportado.nex") #This is to be used in PhyGeo

#2C. STEP - SENSITIVITY PLOT ####

# Read comparisons
comptax_k5  <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k5.csv",  comment.char = "#")
comptax_k10 <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k10.csv", comment.char = "#")
comptax_k15 <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k15.csv", comment.char = "#")
comptax_ew  <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_ew.csv",  comment.char = "#")

# Rename first column to RowID
colnames(comptax_k5)[1] <- "RowID"
colnames(comptax_k10)[1] <- "RowID"
colnames(comptax_k15)[1] <- "RowID"
colnames(comptax_ew)[1] <- "RowID"

# Combine all by RowID
allcomptax <- comptax_k5 %>%
  left_join(comptax_k10, by = "RowID") %>%
  left_join(comptax_k15, by = "RowID") %>%
  left_join(comptax_ew,  by = "RowID")

# Identify the columns that contain "A2"
a2_cols <- grep("A2", colnames(allcomptax), value = TRUE)

# Replace A2 values for "clade_I"
allcomptax[allcomptax$RowID == "clade_I", a2_cols] <-
  allcomptax[allcomptax$RowID == "clade_I_No_Ossturii", a2_cols]

# Replace A2 values for "Osmundaceae"
allcomptax[allcomptax$RowID == "Osmundaceae", a2_cols] <-
  allcomptax[allcomptax$RowID == "Osmundaceae_No_Ossturii", a2_cols]

# Replace A2 values for "Osmundopsis"
allcomptax[allcomptax$RowID == "Osmundopsis", a2_cols] <-
  allcomptax[allcomptax$RowID == "Os_rafaelii_Os_zunigai", a2_cols]

# Discard useless rows 
leftme <- c("Os_rafaelii_Os_zunigai", "Osmundaceae_No_Ossturii", "clade_I_No_Ossturii")
allcomptax <- allcomptax %>%
  filter(!RowID %in% leftme)

#Rename
allcomptax$RowID[1] <- "clade I"
allcomptax$RowID[2] <- "clade leptopteroid"

# Reshape to long format
df_long_allcomptax <- allcomptax %>%
  pivot_longer(
    cols = -RowID,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Category = str_extract(Variable, "k5|k10|k15|ew"),
    Condition = str_extract(Variable, "A\\d+")
  )

# Plot heatmap
sensitivy_plot <- ggplot(df_long_allcomptax, aes(x = Condition, y = fct_rev(RowID), fill = factor(Value))) +
  geom_tile(color = "white") +
  facet_wrap(~Category, ncol = 2) +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  guides(fill = "none") +
  labs(fill = "Value", x = "", y = "", title = "Sensitivity plot") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, face = "italic"),
    strip.text = element_text(size = 8, face = "bold")
  )

#3.  STEP NEW LTT SECTION (see function above)####

#dated_trees <- my_dated_data$dated_trees #used for a posteriori time calibration

dated_trees <- main_tree@phylo #If dated from other software
dated_trees$root.time <- 308.5

my_ltt <- ltt_and_clime(dated_trees)

#Data is in my_ltt$lttclime and my_ltt$clim_ribbon_data
rescale_fun <- gmst_rescaler(my_ltt$lttclime, my_ltt$clim_ribbon_data)

#Plot LTT and paleotemperature
ltt_plot <- ggplot(my_ltt$lttclime, aes(x = time, y = value, color = variable)) +
  geom_rect(data = rect_data, aes(xmin = - xmin, xmax = - xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +
  geom_line(linewidth = 1) + 
  geom_ribbon(data = my_ltt$clim_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
              fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
  geom_ribbon(data = my_ltt$ltt_ribbon_data, aes(x = median.time, ymin = ltt.lower_95, ymax = ltt.upper_95), 
              fill = "#6f7ab7", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
  scale_color_manual(values = c("GMST_50" = "#ee2400",  # Red for GMST_50
                                "GMST_05" = "#ffb09c",  # Light red for GMST_05
                                "GMST_95" = "#ffb09c",  # Light red for GMST_95
                                "median.ltt" = "blue"  # Green for MCMC LTT
  )) + 
  labs(x = "",
       y = "LTT",
       title = "") +
  scale_y_continuous( sec.axis = sec_axis(~ rescale_fun(.), name = "GMST (??C)")) +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12,  face = "bold"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))+
  scale_x_reverse(limits = c(125, 0),
                  breaks = -breaks_vec
  )


# Save as SVG
svglite("ltt_clime.svg", width = 10, height = 6)  # width/height in inches
grid::grid.draw(ggplotGrob(ltt_plot))
dev.off()

#4.  STEP (OPTIONAL) Alternative with ltt phytools ####
indice <- vector(mode = "integer", length = 0)
for (i in 1:length(dated_trees)){
  rf = wRF.dist(my_dated_data$main_tree@phylo, dated_trees[[i]], normalize = T)
  if(rf <= 0.069){
    indice <- append(indice, i, after = length(indice))
  }
}

sampled <- dated_trees[indice]
secltt <- ltt(sampled, plot = F)

# plot the original ltt object
plot(secltt, show.tree=TRUE, lwd=2,
     log.lineages=FALSE, log="y", bty="n", cex.lab=0.9,
     transparency=0.05, axes=FALSE,
     xlab="Ma")

# axes
axis(1, at=h-seq(0,350,by=5), labels=seq(0,350,by=5), las=1, cex.axis=0.8) 
axis(2, las=1, cex.axis=0.8)

# Add rectangles from 55 to 60
rect(xleft = 55, ybottom = par("usr")[3], 
     xright = 60, ytop = par("usr")[4],
     col = rgb(1, 0, 0, alpha = 0.95),  # Semi-transparent red
     border = NA)  # No border


df <- lapply(seq_along(secltt), function(i) {
  data.frame(
    id   = i,
    time = secltt[[i]]$times,
    ltt  = secltt[[i]]$ltt
  )
}) %>%
  bind_rows()


for (i in 1:length(df$time)) {
  df$time[i] <- max(df$time) - df$time[i]
}

# check the result
dim(df)    # should be 28500 x 3
head(df)

ggplot(df, aes(x = time, y = ltt, group = id)) +
  geom_rect(data = rect_data, aes(xmin = - xmin, xmax = - xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +
  theme_minimal() +
  geom_line(alpha = 0.2, colour = "red") +
  labs(x = "Time (Ma)", y = "Number of lineages") +
  scale_x_reverse(breaks= seq(0, 125, by=5))#limits = c(130, 0), breaks = c(125, 66, 56, 33.9, 23, 5.3, 2.6))


#5.  STEP (OPTIONAL) PLOTTING DATED TREE AND SENS_PLOT AND LTT PLOT TOGETHER ####
sensitivity_plot <- sensitivity_plot +
  theme(
    plot.title = element_text(margin = margin(b = 5)),  # reduce bottom margin
    axis.title.x = element_text(margin = margin(t = 5)),  # reduce top margin
    axis.title.y = element_text(margin = margin(r = 5))
  )

ltt_plot <- ltt_plot +
  theme(
    plot.title = element_text(margin = margin(b = 5)),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  )

left_column <- sensitivity_plot / ltt_plot + plot_layout(heights = c(.6,.4))
final_plot <- (left_column | f3_final) + 
  plot_layout(widths = c(0.4, 1.25)); final_plot

ggsave("Figure_2.svg", plot = final_plot, width = 16, height = 8, device = svglite)

#6.  STEP (OPTIONAL) PLOTTING TREES FOR SUPPLEMENTARY MATERIAL ####
#fsize
fsize <- 3.15
#Read all trees
pi5 <- c("pi5_Arbroles1.tre", "pi5_Arbroles2.tre", "pi5_Arbroles3.tre",
         "pi5_Arbroles4.tre", "pi5_Arbroles5.tre", "pi5_Arbroles6.tre")
pi10 <- c("pi10_Arbroles1.tre", "pi10_Arbroles2.tre", "pi10_Arbroles3.tre",
          "pi10_Arbroles4.tre", "pi10_Arbroles5.tre", "pi10_Arbroles6.tre")
pi15 <- c("pi15_Arbroles1.tre", "pi15_Arbroles2.tre", "pi15_Arbroles3.tre",
          "pi15_Arbroles4.tre", "pi15_Arbroles5.tre", "pi15_Arbroles6.tre")
ew <- c("Arbroles1.tre", "Arbroles2.tre", "Arbroles3.tre",
        "Arbroles4.tre", "Arbroles5.tre", "Arbroles6.tre")

pi5_trees  <- vector("list", 6) ; pi5_ggtrees  <- vector("list", 6)
for (i in 1:6) {
  titleme <- paste("A",i,sep = "")
  pi5_trees[[i]] <- ReadTntTree(pi5[[i]])
  if(inherits(pi5_trees[[i]], "multiPhylo")) { pi5_trees[[i]] <- consensus(pi5_trees[[i]]) } 
  pi5_trees[[i]] <- root.phylo(pi5_trees[[i]], outgroup = "Dipteris_conjugata")
  
  # Abbreviate tip labels
  pi5_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi5_trees[[i]]$tip.label)
  
  # Now plot
  pi5_ggtrees[[i]] <- ggtree(pi5_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
    geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
    xlim(0, 20) +
    labs(title = titleme) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
}

upperrow <- pi5_ggtrees[[1]] | pi5_ggtrees[[2]] | pi5_ggtrees[[3]]
lowerrow <- pi5_ggtrees[[4]] | pi5_ggtrees[[5]] | pi5_ggtrees[[6]]
patchtitle <- "k5"
FigS1_plot <- (upperrow / lowerrow)  + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
FigS1_plot
ggsave("FigS1.svg", plot = FigS1_plot, width = 10, height = 15, device = svglite)

pi10_trees <- vector("list", 6); pi10_ggtrees <- vector("list", 6)
for (i in 1:6) {
  titleme <- paste("A",i,sep = "")
  pi10_trees[[i]] <- ReadTntTree(pi10[[i]])
  if(inherits(pi10_trees[[i]], "multiPhylo")) { pi10_trees[[i]] <- consensus(pi10_trees[[i]]) } 
  pi10_trees[[i]] <- root.phylo(pi10_trees[[i]], outgroup = "Dipteris_conjugata")
  
  # Abbreviate tip labels
  pi10_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi10_trees[[i]]$tip.label)
  
  # Now plot
  pi10_ggtrees[[i]] <- ggtree(pi10_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
    geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
    xlim(0, 20) +
    labs(title = titleme) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  
}
upperrow <- pi10_ggtrees[[1]] | pi10_ggtrees[[2]] | pi10_ggtrees[[3]]
lowerrow <- pi10_ggtrees[[4]] | pi10_ggtrees[[5]] | pi10_ggtrees[[6]]
patchtitle <- "k10"
FigS2_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
FigS2_plot
ggsave("FigS2.svg", plot = FigS2_plot, width = 10, height = 15, device = svglite)


pi15_trees <- vector("list", 6); pi15_ggtrees <- vector("list", 6)
for (i in 1:6) {
  titleme <- paste("A",i,sep = "")
  pi15_trees[[i]] <- ReadTntTree(pi15[[i]])
  if(inherits(pi15_trees[[i]], "multiPhylo")) { pi15_trees[[i]] <- consensus(pi15_trees[[i]]) } 
  pi15_trees[[i]] <- root.phylo(pi15_trees[[i]], outgroup = "Dipteris_conjugata")
  
  # Abbreviate tip labels
  pi15_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi15_trees[[i]]$tip.label)
  
  # Now plot
  pi15_ggtrees[[i]] <- ggtree(pi15_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
    geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
    xlim(0, 20) +
    labs(title = titleme) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
}
upperrow <- pi15_ggtrees[[1]] | pi15_ggtrees[[2]] | pi15_ggtrees[[3]]
lowerrow <- pi15_ggtrees[[4]] | pi15_ggtrees[[5]] | pi15_ggtrees[[6]]
patchtitle <- "k15"
FigS3_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
FigS3_plot
ggsave("FigS3.svg", plot = FigS3_plot, width = 10, height = 15, device = svglite)


ew_trees   <- vector("list", 6); ew_ggtrees <- vector("list", 6)
for (i in 1:6) {
  titleme <- paste("A",i,sep = "")
  ew_trees[[i]] <- ReadTntTree(ew[[i]])
  if(inherits(ew_trees[[i]], "multiPhylo")) { ew_trees[[i]] <- consensus(ew_trees[[i]]) } 
  ew_trees[[i]] <- root.phylo(ew_trees[[i]], outgroup = "Dipteris_conjugata")
  
  # Abbreviate tip labels
  ew_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", ew_trees[[i]]$tip.label)
  
  # plot
  ew_ggtrees[[i]] <- ggtree(ew_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
    geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
    xlim(0, 20) +
    labs(title = titleme) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  
}
upperrow <- ew_ggtrees[[1]] | ew_ggtrees[[2]] | ew_ggtrees[[3]]
lowerrow <- ew_ggtrees[[4]] | ew_ggtrees[[5]] | ew_ggtrees[[6]]
patchtitle <- "ew"
FigS4_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
FigS4_plot
ggsave("FigS4.svg", plot = FigS4_plot, width = 10, height = 15, device = svglite)

#7.  STEP (OPTIONAL) ALTERNATIVE PLOTTING OF TWO CONTINUOUS VARIABLES ####
ggtree(main_tree, aes(color=body_mass), continuous = 'colour', yscale = "trait") + 
  scale_color_gradient(name = "Body Mass", 
                       low = "blue", high = "red")

main_tree@data$RGA <- main_tree@data$median_trait 

# option 0: two opposite rectangular trees
p1 <- ggtree(main_tree, aes(color = body_mass), size = 1, ladderize = T, right = T) + 
  scale_color_continuous(low = 'blue', high = 'red', name="Body mass") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 125, by = 25),
                     labels = function(x) 125 - x,  # Convert to Ma
                     limits = c(-1, 155))+
  theme(legend.position = "left",
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())   # legend on the left

p2 <- ggtree(main_tree, aes(color = RGA), size = 1, ladderize = T, right = T) + 
  scale_color_continuous(low = 'red', high = 'blue', name="RGA") +
  scale_x_reverse(breaks = seq(0, 125, by = 25),
                  labels = function(x) 125 - x,  # Same time conversion
                  limits = c(155, -1)) +  # Note: limits are reversed too!
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") # legend on the right

labels <- ggtree(main_tree, alpha=0, ladderize = T, right = T) + geom_tiplab(size=2.95, hjust = 0) + xlim(-1,1)
labels$data$x <- -1

# Place them side by side
p1 + labels + p2 + plot_layout(widths = c(2.4,0.8,2.4))
ggsave(filename = "rgabm.svg", device = "svg", width = 40, height = 25, units = "cm")


# option 1: two opposite rectangular trees
# Left tree (normal orientation)
p1 <- ggtree(main_tree, aes(color = body_mass), size = 1, ladderize = T, right = T) + 
  scale_color_continuous(low = 'blue', high = 'red', name = "Body mass") +
  #scale_color_distiller(palette = "Reds", name = "Body mass")+
  geom_tiplab(align = F, hjust = -0.0025, size = 3, fontface = "italic") +  # Labels on the right
  scale_x_continuous(breaks = seq(0, 125, by = 25),
                     labels = function(x) 125 - x,
                     limits = c(-1, 155)) +
  theme_minimal() +
  theme(legend.position = "left",
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Right tree (reversed orientation) - FIXED
p2 <- ggtree(main_tree, aes(color = RGA), size = 1, ladderize = T, right = T) + 
  scale_color_continuous(low = 'red', high = 'blue', name = "RGA")+
  geom_tiplab(align = F, hjust = 1, size = 3, fontface = "italic") +  # Labels on the LEFT for reversed tree
  scale_x_reverse(breaks = seq(0, 125, by = 25),
                  labels = function(x) 125 - x,
                  limits = c(155, -1)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right")

# Place them side by side
p1 + p2 + plot_layout(widths = c(1, 1))
ggsave(filename = "rgabm_2.svg", device = "svg", width = 40, height = 25, units = "cm")





#option 2
ggtree(main_tree, aes(color = body_mass),
       continuous = 'colour', yscale = "trait", lineend="round", linejoin="round", size=1, alpha=1) + 
  scale_color_gradient(name = "Body Mass", 
                       low = "blue", high = "red") + 
  geom_tiplab(align = T, size = 2.5, colour="black")+
  theme_minimal() + 
  new_scale_color() +
  geom_tree(data = main_tree, aes(color = RGA), 
            continuous = 'colour', 
            position = position_nudge(x = -0.9),
            alpha = 1,
            size = 1.5,
            lineend = "round",
            linejoin = "round") +
  scale_color_distiller(name = "RGA", 
                        palette = "RdPu",
                        direction = 1,
                        guide = guide_colorbar(override.aes = list(alpha = 1))) +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks = seq(0, 125, by = 25),
                     labels = function(x) 125 - x,
                     limits = c(-1, 155))






###* ************************** BLOCK 2 ************************** ####
###* **************************   DTT ***************************   *#
###* ************************************************************  ###

#8.  STEP - DISPARITY THROUGH TIME (USES AN INPUT DATED TREE AND MORPHO MATRIX) ####
ages <- c(255, 205, 193, 65, 49)
input_tree <- "dated_tree_input.tre"

michi <- Claddis::read_nexus_matrix("leptop_disp.nex")
micho <- my_htu_dataset(input_tree = input_tree, input_data = michi, from_file = T)
gas <- process_disp(micho, ages)
cairo_pdf(filename = "DTT-all.pdf", width = 5, height = 7)
plot_my_plate(
  gas,
  title = "Disparity-through-time and phylomorphospace",
  sqformat = FALSE,
  vformat = TRUE)
dev.off()

#8B. STEP - Phylomorphospace
#check node numbers
input_tree <- "dated_tree_input.tre"
tree <- read.mrbayes(input_tree) 
ggtree(tree@phylo, layout = "rectangular", ladderize=TRUE, right=TRUE,
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2, linesize = 0.01, color="black",  offset = 0.5)+
  geom_label(aes(label=node), size=2, color="purple", position = "dodge")

list_A <- unlist(Descendants(tree@phylo, node = 48, type = "tips"))
list_A <- tree@phylo$tip.label[list_A]
list_B <- unlist(Descendants(tree@phylo, node = 63, type = "tips"))
list_B <- tree@phylo$tip.label[list_B]

taxon_groups <- list(
  "osmundeoid" = list_A,
  "leptopteroid" = list_B 
  );class(taxon_groups) <- "taxonGroups"

mat <- read_nexus_matrix("leptop_disp.nex")
tree@phylo$root.time <- 300
matinput <- ordinate_cladistic_matrix(cladistic_matrix = mat,
                                      estimate_all_nodes = T,
                                      time_tree = tree@phylo)

cairo_pdf(filename = "phylomorphospace.pdf", width = 6, height = 5)
plot_morphospace(matinput,
                 taxon_groups = taxon_groups, 
                 plot_group_legend = T,
                 group_legend_position = "bottom_left",
                 z_axis = 3,
                 plot_internal_nodes = F,
                 plot_edges = T,
                 plot_root = T,
                 root_colour = "red", 
                 plot_convex_hulls = T,
                 plot_taxon_names = F,
                 y_limits=c(-1.1,1.25),
                 x_limits = c(-1, 1),
                 palette = "Dark2") 
dev.off()

plot_multi_morphospace(matinput, n_axes = 3,
                       taxon_groups = taxon_groups, 
                       plot_group_legend = T,
                       plot_taxon_names = F,
                       plot_internal_nodes = F,
                       plot_edges = F,
                       plot_root = T,
                       root_colour = "red", 
                       plot_convex_hulls = T,
                       palette = "Dark2");n_axes=2

jk <- plot_chronophylomorphospace(matinput,
                            taxon_groups = taxon_groups)




###* ************************** BLOCK 3 ****************************** ####
###* ********************* EXTRACT ??C AND PLOT REVBAYES ************     *#
###* *****************************************************************  ###

#9. STEP (Optional) - EXTRACT (AND PRINT REVBAYES COMMANDS) °C AND TIME INTERVALS FROM PALEOCLIMATE DATASET (STARTING FROM 315)####
my_limits <- extract_temp(limits = c(325,0))

#10. STEP - PLOTTING OUTPUT FROM MRBAYES SFBD #####
## Import summary tree produced by Mr. Bayes
input_tree <- "dated_tree_input.tre"

tree <- read.mrbayes(input_tree) 

scaff <- ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
                size = 2.05)+
  geom_tiplab(size=3, color="black",  offset = 0.1)


scaff + geom_label(aes(label=node), size=2, color="purple", position = "dodge")
scaff

## Import all log (.p) files using EvoPhylo
Comb_posterior <- combine_log(".", burnin = 0.25, downsample = 2500)

#Extract speciation and extinction columns
pertime_rates <- Comb_posterior[,11:26]
pertime_rate_speciation <- pertime_rates[,1:8]
pertime_rate_extinction <- pertime_rates[,9:ncol(pertime_rates)]

#Reshape to long format
pertime_long_sp <- pertime_rate_speciation %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "time_period"),
    names_pattern = "(.*)_(\\d+)"
  )
colnames(pertime_long_sp)[2]  <- "value" 
pertime_long_sp$type <- "net speciation"

#Reshape to long format
pertime_long_ext <- pertime_rate_extinction %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "time_period"),
    names_pattern = "(.*)_(\\d+)"
  )
colnames(pertime_long_ext)[2]  <- "value" 
pertime_long_ext$type <- "relative extinction"

pertime_long <- rbind(pertime_long_sp, pertime_long_ext)

pertime_long <- pertime_long %>%
  mutate(
    time_period = recode(
      time_period,
      "1" = "TBi: Root-193",
      "2" = "TBii: 193-49",
      "3" = "TBiii: 49-16",
      "4" = "TBiv: 16-0")
  )


#Plots in the following point
ratesthroughtime <- ggplot(pertime_long, aes(x = value, fill = type)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 25, 
                 color = "black") +
  facet_wrap(~ time_period, scales = "free") +
  #facet_wrap(~ facet_label, scales = "free")+
  labs(
    title = "Distribution of Evolutionary Rates Across Time Periods",
    x = "Rate Value",
    y = "Frequency",
    fill = "Rate type"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgray"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )


#10B. Fig S5 in Urrea et al 2026 ####
input_tree <- "dated_tree_input.tre"
tree <- read.mrbayes(input_tree)

scaff <- ggtree(
  tree,
  layout = "rectangular",
  ladderize = TRUE,
  right = TRUE
)


## Node labels to annotate
labels_df <- tibble(
  node  = c(49, 55, 63, 68, 72, 69),
  label = c("12", "25", "42", "59", "62", "65")
)

node_labels <- scaff$data %>%
  mutate(node = as.integer(node)) %>%
  select(node, x, y) %>%
  inner_join(labels_df, by = "node")

## 10. Final plot
scaff +
  geom_tree(size = 1.25) +
  geom_tiplab(size = 3, colour = "black", offset = 0.1) +
  geom_text(data = node_labels, aes(x = x, y = y, label = label),
    size = 2.75, hjust = 1.25, vjust = -0.25, fontface = "bold", colour = "black") +
  scale_x_continuous(breaks = breaks_vec, labels = abs(breaks_vec), expand = expansion(mult = c(0.05, 0.38)))
ggsave(filename = "Fig_S5A.svg", device = "svg", height = 20, width = 15, units = "cm")
ggsave(filename = "Fig_S5A.png", device = "png", height = 20, width = 15, units = "cm")


#11. STEP - Rates and selection type (mostly edited from Simoes et al 2022 paper) ####
#Input tree
input_tree <- "dated_tree_input.tre"
tree <- read.mrbayes(input_tree)
dated_trees <- tree@phylo


#LTT
my_ltt <- ltt_and_clime(dated_trees) #from step 3 function

#Data is in my_ltt$lttclime and my_ltt$clim_ribbon_data
rescale_fun <- gmst_rescaler(my_ltt$lttclime, my_ltt$clim_ribbon_data)


#Plot LTT and paleotemperature
ltt_plot <- ggplot(my_ltt$lttclime, aes(x = time, y = value, color = variable)) +
  geom_rect(data = rect_data, aes(xmin = - xmin, xmax = - xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +
  geom_line(linewidth = 1) + 
  geom_ribbon(data = my_ltt$clim_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
              fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
  scale_color_manual(values = c("GMST_50" = "#ee2400",  # Red for GMST_50
                                "GMST_05" = "#ffb09c",  # Light red for GMST_05
                                "GMST_95" = "#ffb09c",  # Light red for GMST_95
                                "median.ltt" = "blue"  # Green for MCMC LTT
  )) + 
  labs(x = "Ma",
       y = "LTT",
       title = "Lineages through time and global mean surface temperature") +
  scale_y_continuous( sec.axis = sec_axis(~ rescale_fun(.), name = "GMST")) +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))+
  scale_x_reverse(limits = c((breaks_vec[1]*(-1)), 0.05),
                  breaks = -breaks_vec
  )
ltt_plot




#check node numbers
ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2, linesize = 0.01, color="black",  offset = 0.5)+
  geom_label(aes(label=node), size=2, color="purple", position = "dodge")


## Import all log (.p) files using EvoPhylo
Comb_posterior <- combine_log(".", burnin = 0.25, downsample = 2500)

## Get table of clock rates with summary stats for each node in the tree for each relaxed clock partition
RateTable_Medians_no_clades <- get_clockrate_table_MrBayes(tree, summary = "median")
node_vec <- RateTable_Medians_no_clades$nodes

#Check the actual number in the plot
Glecheniales <- Descendants(tree@phylo,
                            node = 85,
                            type = "all")
Glecheniales <- append(Glecheniales, 85, after = length(Glecheniales))

Osmunda <- Descendants(tree@phylo,
                       node = 49,
                       type = "all")
Osmunda <- append(Osmunda, 49, after = length(Osmunda))

Osmundopsis <- Descendants(tree@phylo,
                           node = 55,
                           type = "all")
Osmundopsis <- append(Osmundopsis, 55, after = length(Osmundopsis))

Leptoperoid <- Descendants(tree@phylo, node = 63, type = "all")
Leptoperoid <- append(Leptoperoid, 63, after = length(Leptoperoid))


node_vec[node_vec %in% Osmunda] <- "Osmunda"
node_vec[node_vec %in% Glecheniales] <- "Glecheniales"
node_vec[node_vec %in% Leptoperoid] <- "Leptoperoid"
node_vec[node_vec %in% Osmundopsis] <- "Osmundopsis"
node_vec[node_vec %in% 48] <- "Osmundopsis+Osmunda"
node_vec[node_vec %in% 46] <- "Osmundaceae"
#node_vec[node_vec %in% 34] <- "T.muelleri"
node_vec[node_vec %in% 47] <- "Osmundaceae[-Osmundastrum]"
#node_vec[node_vec %in% 42] <- "crown Osmundaceae"
node_vec[node_vec %in% 45] <- "Polypodiopsida"
node_vec[node_vec %in% 5] <- "Osmundastrum"

clade <- node_vec

cladeMedians <- cbind(clade,RateTable_Medians_no_clades)

## Import rate table with clade membership (new "clade" column added)
RateTable_Medians1 <- cladeMedians

## Get summary statistics table for each clade by clock
clockrate_summary(RateTable_Medians1, digits=2)

## Transform table from wide to long format
RatesByClade <- clock_reshape(RateTable_Medians1)

## Get table of pairwise t-tests for difference between the posterior mean and the rate for each tree node
get_pwt_rates_MrBayes(RateTable_Medians1, Comb_posterior)

Sel <- plot_treerates_sgn(
  type = "MrBayes", tree, Comb_posterior,
  trans = "log",
  clock = 1,
  summary = "mean",
  branch_size = 1, tip_size = 2.2,
  xlim = c(-375, (Ntip(tree) + 25)),
  nbreaks = 10,
  geo_size = list(2, 2),
  geo_skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
  threshold = c("90%", "95%"),
  low = "blue", mid = "purple", high = "red"
)

object <- Sel$data

# Leptopteroid ancestor
lepclade <- getMRCA(
  phy = as.phylo(object),
  tip = c("Cacumen_expansa", "Todites_muelleri")
)

# Base ggtree object and coordinates
base_tree <- ggtree(object, ladderize = TRUE, right = TRUE)
df_coords <- base_tree$data

# Make sure node is integer
df_coords$node <- as.integer(df_coords$node)

# node labels
labels_df <- data.frame(
  node  = c(46, 48, 49, 55, 63, 82, 68, 74),
  label = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII"),
  stringsAsFactors = FALSE
)

# Join labels with coordinates
node_labels <- df_coords %>%
  select(node, x, y) %>%
  inner_join(labels_df, by = "node")

# Actual plot
A0 <- base_tree +
  theme_tree2() +  
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +
  geom_tree(layout = "rectangular",
    size = 1.25, colour = "#5E5D5D") +
  geom_tiplab(size = 3.0, offset = 1.75, fontface = "bold.italic") +
  geom_text(data = node_labels, aes(x = x, y = y, label = label), fontface = "bold", color = "darkgrey",
    size = 2.75, hjust = 2.5, vjust = -1.15) +
  geom_text(aes(x = branch, label = ifelse(round(age_median, 0) > 0 & !isTip, round(age_median, 0), "")),
    size = 2.75, hjust = 0.8, vjust = 1.75) +
  geom_text(aes(x = branch, label = ifelse(round(prob, 1) > 0.5 & !isTip, round(prob, 1), "")),
    size = 2.75, hjust = 1, vjust = -1.05, fontface = "italic") +
  theme(legend.position = "bottom",
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9)) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 0))) +
  scale_color_brewer(palette = "Dark2") +
  geom_cladelab(
    node = lepclade,
    label = "leptopteroids",
    align = FALSE,
    offset = 140,
    offset.text = 6,
    angle = -90,
    fontsize = 4
  ) +
  geom_range(
    range = "age_0.95HPD",
    colour = "red",
    size = 1.25,
    alpha = 0.25,
    position = position_nudge(x = -0.5)
  ) +
  coord_geo(
    xlim = c(breaks_vec[1], Ntip(tree) + 1),
    ylim = c(0.1, max(object$y)),
    expand = TRUE,
    dat = list("eras", "periods", "epochs"),
    abbrv = list(FALSE, FALSE, TRUE),
    skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
    pos = list("top", "top", "top"),
    alpha = 1,
    height = unit(0.75, "line"),
    size = "auto",
    fittext_args = list(size = 22),
    rot = 0,
    neg = TRUE
  ) +
  scale_x_continuous(
    breaks = breaks_vec,
    labels = abs(breaks_vec),
    expand = expansion(mult = c(0.05, 0.38))
  )
A0

#Fig. 1
panelB <- ratesthroughtime / ltt_plot #ratesthroughtime is taken in from previous steps
A0 | panelB + plot_layout(widths = c(2, 1))  # A0 gets 2/3, panelB gets 1/3
ggsave("Fig1.svg", device = "svg", width = 37, height = 27, units = "cm")
ggsave("Fig1.eps", device = "eps", width = 37, height = 27, units = "cm")



#NOW with the selecction mapped ###
# Mapping from numeric clockfac descriptive category
selection_labels <- c(
  "1" = "stabilising (below 95 CI)",
  "2" = "stabilising-neutral (lower 90???95 CI bound)",
  "3" = "neutral (90 CI)",
  "4" = "positive-neutral (upper 90-95 CI bound)",
  "5" = "positive (above 95 CI)"
)

# Fixed, consistent, colourblind-safe palette
selection_cols <- c(
  "stabilising (below 95 CI)"          = "#0978d9",
  "stabilising-neutral (lower 90???95 CI bound)"     = "#2e948d",
  "neutral (90 CI)"                    = "#96ab90",
  "positive-neutral (upper 90-95 CI bound)"        = "#d9642e",
  "positive (above 95 CI)"             = "#b0160e"
)

#Partition 1
Sel <- plot_treerates_sgn(
  type = "MrBayes", tree, Comb_posterior,
  trans = "log",
  clock = 1,
  summary = "mean",
  branch_size = 1, tip_size = 2.2,
  xlim = c(-375, (Ntip(tree) + 25)),
  nbreaks = 10,
  geo_size = list(2, 2),
  geo_skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
  threshold = c("90%", "95%"),
  low = "blue", mid = "purple", high = "red"
)

object <- Sel$data

# Attach new factor with fixed levels
object$selection <- factor(
  selection_labels[as.character(object$clockfac)],
  levels = selection_labels)

A1 <- ggtree(object, ladderize = TRUE, right = TRUE) +
  theme_tree2() + 
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +  
  geom_tree(layout = "rectangular", size = 1.5,
            aes(colour = selection)) +  
  geom_tiplab(size = 3, offset = 1.75,
              fontface = "bold.italic", aes(colour = selection)) +  
  theme(legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9)) +  
  guides(colour = guide_legend(
    override.aes = list(shape = 16, size = 4)
  )) +  
  scale_color_manual(values = selection_cols, drop = FALSE) +
  coord_geo(
    xlim = c(inferior_limit, Ntip(tree) + 1),
    ylim = c(0.1, max(object$y)), expand = TRUE,
    dat = list("eras", "periods", "epochs"),
    abbrv = list(FALSE, FALSE, TRUE),
    skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
    pos = list("top", "top", "top"),
    alpha = 1, height = unit(0.75, "line"),
    size = "auto", fittext_args = list(size = 22),
    rot = 0, neg = TRUE
  ) +
  scale_x_continuous(
    breaks = breaks_vec,
    labels = abs(breaks_vec),
    expand = expansion(mult = c(0.001, 0.15))
  )
A1
ggsave(filename = "selection_vegetative.svg", device = "svg", width = 20, height = 15, units = "cm")


#Partition 2
Sel <- plot_treerates_sgn(
  type = "MrBayes", tree, Comb_posterior,
  trans = "log",
  clock = 2,
  summary = "mean",
  branch_size = 1, tip_size = 2.2,
  xlim = c(-375, (Ntip(tree) + 25)),
  nbreaks = 10,
  geo_size = list(2, 2),
  geo_skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
  threshold = c("90%", "95%"),
  low = "blue", mid = "purple", high = "red"
)

object <- Sel$data

object$selection <- factor(
  selection_labels[as.character(object$clockfac)],
  levels = selection_labels
)

A2 <- ggtree(object, ladderize = TRUE, right = TRUE) +
  theme_tree2() +  
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rgbfill, alpha = my_alpha, inherit.aes = FALSE) +  
  geom_tree(layout = "rectangular", size = 1.5,
            aes(colour = selection)) +
  geom_tiplab(size = 3.0, offset = 1.75,
              fontface = "bold.italic",
              aes(colour = selection)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9)) +  
  guides(colour = guide_legend(
    override.aes = list(shape = 16, size = 4)
  )) +  
  scale_color_manual(values = selection_cols, drop = FALSE) +  
  coord_geo(
    xlim = c(inferior_limit, Ntip(tree) + 1),
    ylim = c(0.1, max(object$y)), expand = TRUE,
    dat = list("eras", "periods", "epochs"),
    abbrv = list(FALSE, FALSE, TRUE),
    skip = c("Pliocene", "Pleistocene", "Holocene", "Quaternary"),
    pos = list("top", "top", "top"),
    alpha = 1, height = unit(0.75, "line"),
    size = "auto", fittext_args = list(size = 22),
    rot = 0, neg = TRUE
  ) +  
  scale_x_continuous(
    breaks = breaks_vec,
    labels = abs(breaks_vec),
    expand = expansion(mult = c(0.001, 0.15))
  )
A2
ggsave(filename = "selection_fertile.svg", device = "svg", width = 20, height = 15, units = "cm")


#cowplot::plot_grid(A1, A2, ncol = 2, rel_widths = c(1, 1))


#DTT FOR SUPPLEMENTARY MATERIAL ####
vegchar <- c(1, 2, 3, 5, 6, 7, 14, 15, 16, 17, 18, 19, 20, 21, 26, 31) #This the list of vegetative characters
fertchar <- c(4, 8, 9, 10, 11, 12, 13, 22, 23, 24, 25, 27, 28, 29, 30) #This is the list of fertile characters
mytreeinput <- "dated_tree_input.tre"

veg <- splitme(mydata = michi, bytype = F, chars = vegchar)$custom #splitting matrix by vegs
veghtu <- my_htu_dataset(input_tree = mytreeinput, input_data = veg, from_file = T)
vegdisp <- process_disp(veghtu)

fer <- splitme(mydata = michi, bytype = F, chars = vegchar)$custom #splitting matrix by vegs
ferhtu <- my_htu_dataset(input_tree = mytreeinput, input_data = fer, from_file = T)
ferdisp <- process_disp(ferhtu)

nominal <- c(1, 6, 7, 16, 25, 27, 28) 
ordinal <- c(2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 26, 29, 30, 31)
nominal_data <- splitme(mydata = michi, bytype = F, chars = nominal)$custom #splitting matrix by vegs
nominalhtu <- my_htu_dataset(input_tree = mytreeinput, input_data = nominal_data, from_file = T)
nominaldisp <- process_disp(nominalhtu)

ordinal_data <- splitme(mydata = michi, bytype = F, chars = ordinal)$custom #splitting matrix by vegs
ordinalhtu <- my_htu_dataset(input_tree = mytreeinput, input_data = ordinal_data, from_file = T)
ordinaldisp <- process_disp(ordinalhtu)

# These generated too much NA 
# nominal_veg <- c(0, 5, 6, 15)
# nominal_fer <- c(24, 26, 27)
# nominalveg_data <- splitme(mydata = michi, bytype = F, chars = nominal_veg)$custom #splitting matrix by vegs
# nominalveghtu <- my_htu_dataset(input_tree = mytreeinput, input_data = nominalveg_data, from_file = T)
# nominalvegdisp <- process_disp(nominalveghtu)
# nominalfer_data <- splitme(mydata = michi, bytype = F, chars = nominal_fer)$custom #splitting matrix by vegs
# nominalferhtu <- my_htu_dataset(input_tree = mytreeinput, input_data = nominalfer_data, from_file = T)
# nominalferdisp <- process_disp(nominalferhtu)

ordinal_veg <- c(1, 2, 3, 4, 13, 14, 16, 17, 18, 19, 20, 25, 30)
ordinal_fer <- c(7, 8, 9, 10, 11, 12, 21, 22, 23, 28, 29)
ordinalveg_data <- splitme(mydata = michi, bytype = F, chars = ordinal_veg)$custom #splitting matrix by vegs
ordinalveghtu <- my_htu_dataset(input_tree = mytreeinput, input_data = ordinalveg_data, from_file = T)
ordinalvegdisp <- process_disp(ordinalveghtu)

ordinalfer_data <- splitme(mydata = michi, bytype = F, chars = ordinal_fer)$custom #splitting matrix by vegs
ordinalferhtu <- my_htu_dataset(input_tree = mytreeinput, input_data = ordinalfer_data, from_file = T)
ordinalferdisp <- process_disp(ordinalferhtu)


plot_my_plate(vegdisp, title="Vegetative partition", vformat=T, sqformat=F)
plot_my_plate(ferdisp, title="Fertile partition", vformat=T, sqformat=F)
plot_my_plate(nominaldisp, title="Nominal partition", vformat=T, sqformat=F)
plot_my_plate(ordinaldisp, title="Ordinal partition", vformat=T, sqformat=F)
plot_my_plate(ordinalvegdisp, title="Ordinal vegetative partition", vformat=T, sqformat=F)
plot_my_plate(ordinalferdisp, title="Ordinal fertile partition", vformat=T, sqformat=F)



#TREES FROM ALL MODELS FOR SUPP. INFO####
statedpartitioned_trees <- c("leptop_igrdiv.nex.con.tre",
                             "leptop_igrran.nex.con.tre",
                             "leptop_tk02ran.nex.con.tre",
                             "leptop_tk02div.nex.con.tre")

plot_multiple_mrbayes_trees(statedpartitioned_trees,
                            ncol = 2,
                            titles = c("IGR Div", "IGR Ran", "TK02 Ran", "TK02 Div"),
                            filename = "combined_mrbayes_trees_A.svg")


sourcepartitioned_trees <- c("leptop_igrran_orgpart.nex.con.tre",
                             "leptop_igrdiv_orgpart.nex.con.tre",
                             "leptop_tk02ran_orgpart.nex.con.tre",
                             "leptop_tk02div_orgpart.nex.con.tre")

plot_multiple_mrbayes_trees(sourcepartitioned_trees,
                            ncol = 2,
                            titles = c("IGR Ran by source", "IGR Div by source", "TK02 Ran by source", "TK02 Div by source"),
                            filename = "combined_mrbayes_trees_B.svg")

fossiltip_trees <- c("leptop_tk02fossiltip.nex.con.tre",
                     "leptop_tk02fossiltip.nex.con.tre",
                     "leptop_igrfossiltip.nex.con.tre",
                     "leptop_igrfossiltip.nex.con.tre")

plot_multiple_mrbayes_trees(fossiltip_trees,
                            ncol = 2,
                            titles = c("Fossiltip TK02 by state", "Fossiltip TK02 by source", "Fossiltip IGR by state", "Fossiltip IGR by source"),
                            filename = "fossiltip_trees.svg")


iw5 <- ReadTntTree("iw5.tre"); write.beast(iw5, file = "iw5.nex", translate = F)
iw10 <- ReadTntTree("iw10.tre"); write.beast(iw10, file = "iw10.nex", translate = F)
iw15 <- ReadTntTree("iw15.tre"); write.beast(iw15, file = "iw15.nex", translate = F)

fossiltip_trees <- c("iw5.nex",
                     "iw10.nex",
                     "iw15.nex",
                     "leptop_mb.nex.con.tre")

plot_multiple_mrbayes_trees(fossiltip_trees,
                            ncol = 2,
                            titles = c("Implied weighting K5",
                                       "Implied weighting K10",
                                       "Implied weighting K15",
                                       "Non-clock BI"),
                            filename = "nonclock_trees.svg", nexfile = TRUE)
