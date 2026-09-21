# This function evaluates a MaxEnt species distribution model using k-fold cross-validation. 
# Occurrence records are divided into training and test datasets according to predefined fold assignments. 
# For each fold, a MaxEnt model is fitted using the training occurrences and background points, 
# while the test occurrences are used for model evaluation. 

# Arguments:
# p                - Data frame containing occurrence coordinates in the first two columns and cross-validation fold assignments in the third column. 
# env              - RasterStack or RasterBrick containing environmental predictors (from 'raster' package). 
# bg               - Background points used for MaxEnt model fitting. 
# removeDuplicates - Logical value indicating whether duplicate occurrence records should be removed by MaxEnt. 
# maxent.args      - Additional arguments passed to the MaxEnt algorithm. 
# print.progress   - Logical value indicating whether progress information should be printed during cross-validation.

cv.maxent <- function(
    p,
    env,
    bg,
    removeDuplicates = TRUE,
    maxent.args = c(
      "addsamplestobackground=true",
      "autofeature=true",
      "betamultiplier=1"
    ),
    print.progress = TRUE
) {
  # Load required packages
  require(raster)
  require(dismo)
  require(stringi)
  
  # Create a main temporary directory for MaxEnt test data if it doesn't exist
  if (!dir.exists("./tempMaxent")) {
    dir.create("./tempMaxent")
  }
  
  # Create a unique subfolder for this specific run to prevent file conflicts
  temp.folder.name <- stri_rand_strings(1, 10)
  temp_dir_path <- paste0("./tempMaxent/", temp.folder.name)
  dir.create(temp_dir_path)
  
  # Initialize vectors to store training and test AUC values for each fold
  train.auc <- vector("numeric")
  test.auc <- vector("numeric")
  
  # List to store permutation importance for each fold
  perm.imp.list <- list()
    
  # Extract cross-validation fold assignments (assumed to be in the 3rd column)
  folds <- p[, 3]
  max_folds <- max(folds)
  
  # Run k-fold cross-validation
  for (i in seq(1, max_folds)) {
    
    # Print progress if requested
    if (print.progress) {
      print(paste("Cross-validation: Fold", i, "out of", max_folds))
    }
    
    ## 1. Split occurrence data into training and test sets
    train.p <- as.data.frame(p[folds != i, 1:2])
    test.p  <- as.data.frame(p[folds == i, 1:2])
    
    ## 2. Extract environmental variables for test occurrences
    test.swd.env <- as.data.frame(raster::extract(env, test.p))
    
    ## 3. Create a species-with-data (SWD) dataframe for test occurrences
    test.swd <- data.frame(
      Species = "species",
      test.p,
      test.swd.env
    )
    
    # Rename coordinate columns for MaxEnt compatibility
    colnames(test.swd)[2:3] <- c("longitude", "latitude")
    
    # Remove test occurrences that have missing environmental data (NA values)
    test.swd <- test.swd[complete.cases(test.swd), ]
    
    ## 4. Save test occurrences as an SWD file for MaxEnt evaluation
    test_file_path <- paste0(temp_dir_path, "/TestSWD.csv")
    
    write.table(
      test.swd,
      file = test_file_path,
      sep = ",",
      dec = ".",
      row.names = FALSE,
      quote = FALSE
    )
    
    ## 5. Add the test sample file to MaxEnt arguments for this specific fold
    # (Important: use a temporary variable so we don't permanently alter maxent.args for the next iterations)
    current_maxent_args <- c(
      maxent.args,
      paste0("testsamplesfile=", getwd(), substring(test_file_path, 2)) # substring removes the dot from "./"
    )
    
    ## 6. Fit the MaxEnt model using the training occurrences
    mod <- dismo::maxent(
      x = env,
      p = train.p,
      a = bg,
      removeDuplicates = removeDuplicates,
      args = current_maxent_args
    )
    
    ## 7. Extract training and test AUC values from model results
    train.auc[i] <- mod@results["Training.AUC", 1]
    test.auc[i]  <- mod@results["Test.AUC", 1]
    
    # Extract permutation importance values
    res <- mod@results
    
    # Find rows corresponding to permutation importance
    perm_idx <- grep("\\.permutation\\.importance", rownames(res))
    perm_vals <- res[perm_idx, 1]
    
    # Clean variable names
    names(perm_vals) <- gsub("\\.permutation\\.importance", "", names(perm_vals))
    
    # Store in the list
    perm.imp.list[[i]] <- perm_vals
    
    # Clean up the temporary test file for the current fold
    file.remove(test_file_path)
  }
    
  
  # Remove the temporary directory completely after cross-validation finishes
  unlink(temp_dir_path, recursive = TRUE, force = TRUE)
  
  #Calculate mean and sd of permutation importance across all folds
  perm.imp.df <- do.call(rbind, perm.imp.list)
  mean.perm.imp <- colMeans(perm.imp.df, na.rm = TRUE)
  sd.perm.imp <- apply(perm.imp.df, 2, sd, na.rm = TRUE)
  
  # 8. Return cross-validation performance statistics
  return(
    list(
      mean.train.auc   = mean(train.auc, na.rm = TRUE),
      sd.train.auc     = sd(train.auc, na.rm = TRUE),
      
      mean.test.auc    = mean(test.auc, na.rm = TRUE),
      sd.test.auc      = sd(test.auc, na.rm = TRUE),
      
      # Overfitting is calculated as the difference between Training and Test AUC
      mean.overfitting = mean(train.auc - test.auc, na.rm = TRUE),
      sd.overfitting   = sd(train.auc - test.auc, na.rm = TRUE),
      mean.permutation.importance = mean.perm.imp,
      sd.permutation.importance   = sd.perm.imp
    )
  )
}