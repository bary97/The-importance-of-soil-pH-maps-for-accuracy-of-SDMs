# Environmental filtering of species occurrences (Varela et al., 2014) 
# and selection of background points restricted to unoccupied environmental space


# This function creates a multidimensional environmental grid using scaled 
# predictors (with a defined bin width in standard deviations) and performs 
# two main tasks:
# 1. Filters species occurrence records by retaining a single random record 
#    per environmental cell to mitigate sampling bias.
# 2. Optionally filters background points by removing any points 
#    that fall into environmental cells occupied by species occurrences, 
#    followed by random subsampling to the required sample size.


# Arguments:
# occ_data      - Data frame containing species occurrence records with coordinates and extracted scaled environmental variables.
# pred_cols     - Vector of column names for predictors. 
# bg_candidates - Optional data frame of candidate background points (including coordinates and extracted scaled predictors).
# n_bg          - Required number of final background points after filtering.
# step          - Environmental grid cell width (bin size) in standard deviations. Default is 0.1.
# seed          - Random seed for reproducibility. Default is 700.

enviFilter <- function(occ_data, 
                       pred_cols,            
                       bg_candidates = NULL, 
                       n_bg = NULL, 
                       step = 0.1, 
                       seed = 700) {
  
  # Extract values of predictors for occurrences and check NA values
  pred_occ <- as.data.frame(occ_data[, pred_cols, drop = FALSE])
  
  if (anyNA(pred_occ)) {
    stop("Occurrence predictors contain NA values.")
  }
  
  # Create a multidimensional environmental grid based on occurrence ranges
  grid <- lapply(pred_occ, function(x) {
    # floor and ceiling ensure min/max bounds are neatly snapped to exact multiples of the step size  
    min_x <- floor(min(x) / step) * step
    max_x <- ceiling(max(x) / step) * step
    seq(min_x, max_x, by = step)
  })
  
  # Assign occurrence records to discrete environmental grid cells
  occ_groups <- data.frame(ID = seq_len(nrow(occ_data)))
  
  for (i in seq_along(pred_occ)) {
    occ_groups[[paste0("env_", i)]] <- cut(
      pred_occ[[i]],
      breaks = grid[[i]],
      include.lowest = TRUE,
      labels = FALSE
    )
  }
  
  # Create a unique multidimensional cell identifier for each occurrence by combining its bin indices across all predictors
  occ_groups$cell <- apply(occ_groups[, -1, drop = FALSE], 1, paste, collapse = "_")
  
  # Randomly shuffle occurrence records to ensure unbiased selection, retain a single record per unique environmental cell, and extract the filtered subset
  set.seed(seed)
  filtered_occ <- occ_groups[sample(seq_len(nrow(occ_groups))), ]
  filtered_occ <- filtered_occ[!duplicated(filtered_occ$cell), ]
  
  result <- list(occ_filtered = occ_data[filtered_occ$ID, ])
  
  # Filter candidate background points automatically using the same pred_cols
  if (!is.null(bg_candidates)) {
    # Extract values of predictors for background points and check NA values
    pred_bg <- as.data.frame(bg_candidates[, pred_cols, drop = FALSE])
    
    if (anyNA(pred_bg)) {
      stop("Background predictors contain NA values.")
    }
    
    # Identify unique environmental cells occupied by species occurrences and initialise background groups tracking table
    occupied_cells <- unique(occ_groups$cell)
    bg_groups <- data.frame(ID = seq_len(nrow(bg_candidates)))
    
    # Assign background candidates to the discrete environmental grid cells
    for (i in seq_along(pred_bg)) {
      bg_groups[[paste0("env_", i)]] <- cut(
        pred_bg[[i]],
        breaks = grid[[i]],
        include.lowest = TRUE,
        labels = FALSE
      )
    }
    # Create a unique multidimensional cell identifier for each background candidate
    bg_groups$cell <- apply(bg_groups[, -1, drop = FALSE], 1, paste, collapse = "_")
    # Remove background points falling within environmental cells occupied by species occurrences
    bg_groups <- bg_groups[!(bg_groups$cell %in% occupied_cells), , drop = FALSE]
    
    # Verify that an adequate number of background points remains after environmental filtering
    if (!is.null(n_bg) && nrow(bg_groups) < n_bg) {
      stop(paste0("Only ", nrow(bg_groups), " background points remain after environmental filtering, ",
                  "but ", n_bg, " are required. Please increase the number of candidate background points."))
    }
    
    # Randomly sample the required number of background points from the filtered pool without replacement
    set.seed(seed)
    selected_bg_ids <- sample(
      bg_groups$ID,
      size = if (!is.null(n_bg)) n_bg else nrow(bg_groups),
      replace = FALSE
    )
    # Extract the final filtered background records and compile them into the output list
    result$bg_filtered <- bg_candidates[sort(selected_bg_ids), ]
  }
  
  return(result)
}

