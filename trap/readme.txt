TRAP prediction with spatial random forest (spatRF)

References:
code for preprocessing mostly came from Magali's work (https://github.com/magali17/TRAP)

spatRF code was largely based on the R package by Travis (https://github.com/theewai/spatRF)

1. preprocessing for geocovariates
preprocess.R (relies on functions defined by preprocessing_function.R)

2. spatRF
spatTree.R (original work by Travis, relies on aux_functions.R)
spatTreeModified.R (modified by Si, allows for tuning multiple parameters)

3. TRAP models
model.R -- cross-validation to compare the performance of different models
pred_grid.R -- make predictions on grids
full_summary.Rmd -- details on the model, including variable importance measures and visualization