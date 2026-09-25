# The importance of soil pH maps for accuracy of plant distribution models

This repository contains the R scripts used for data processing and species distribution modelling in the paper *"The importance of soil pH maps for accuracy of plant distribution models"*. 

## Repository Contents

The repository includes custom R functions designed to mitigate spatial sampling bias and properly evaluate MaxEnt models using cross-validation:

*   **`enviFilter.R`**: A function for environmental filtering of species occurrences and background points (following Varela et al., 2014). It creates a multidimensional environmental grid to retain a single occurrence record per environmental cell, mitigating sampling bias. It also filters background points by removing those that fall into environmental cells occupied by species occurrences.
*   **`cv.maxent.R`**: A function to evaluate MaxEnt species distribution models using k-fold cross-validation. It calculates model performance metrics including mean and standard deviation of training AUC, test AUC, model overfitting (difference between training and test AUC), and predictor permutation importance.

## Dependencies

The scripts require the following R packages:
*   `raster`
*   `dismo`
*   `stringi`

## Data Availability

The spatial data and initial occurrence dataset used by these scripts are too large to be hosted on GitHub. They are permanently archived and freely available on Zenodo.

To run the analyses, download the data from Zenodo and place the files into a local `/data` directory before executing the scripts.
