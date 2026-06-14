Random Forest classification of successional stages in subtropical Atlantic Forest

This repository contains the R code supporting the study "Random Forest classification of successional stages in subtropical Atlantic Forest: integrating stand structure and floristic composition", which develops a quantitative, data-driven framework to classify successional stages (early, medium, advanced) in the Southern Brazilian Atlantic Forest of Santa Catarina.

A Random Forest classifier is trained and validated on systematic forest-inventory data, integrating stand structural attributes with floristic-composition gradients derived from a Principal Coordinates Analysis (PCoA). The resulting classifications are then compared against the legal structural thresholds of CONAMA Resolution 04/94.

Pipeline overview

The script random_forest_succession.R runs the full analytical workflow in sequence:


Data import — loads the tree, metadata and stand-parameter tables and merges them by plot.
Pre-processing — inspects and imputes missing values (median imputation).
Collinearity removal (VIF) — iterative, stepwise pruning of predictors by Variance Inflation Factor and pairwise Pearson correlation.
Variable selection (Boruta) — confirms the informative predictors against shadow features.
Floristic gradients (PCoA) — ordination of a Bray–Curtis dissimilarity matrix; the first ten axes are retained as compositional predictors.
Train/test split — stratified 70/30 partition by successional stage.
ntree optimisation — selection of the number of trees from the Out-Of-Bag error curve.
Random Forest model — training of the final classifier.
Evaluation — confusion matrices, per-class metrics, ROC curves and multi-class AUC.
Model interpretation — Partial Dependence Plots for the most important predictors.
Comparison with CONAMA 04/94 — observed and predicted structure contrasted against the legal thresholds, with a compliance assessment.


Requirements

The analysis was developed in R. The following packages are required:

rinstall.packages(c(
  "readxl", "dplyr", "tidyr", "vegan", "ape", "randomForest",
  "caret", "ggplot2", "Boruta", "pROC", "corrplot", "car",
  "janitor", "pdp", "gridExtra", "cowplot"
))

Usage

Place the input data files inside a data/ folder in the project root, then run:

rsource("random_forest_succession.R")

The script prints progress and summary statistics to the console and writes the figures (TIFF, 300 dpi) to the working directory.

Data availability

The raw data of the Santa Catarina Forest Inventory (IFFSC) are not included in this repository. They may be requested from the State of Santa Catarina through the official channels of the inventory (https://iff.sc.gov.br). The code is provided so that the analytical workflow can be inspected and reproduced once access to the data has been granted.

Citation

If you use this code, please cite the associated article. Full citation details will be added upon publication.
