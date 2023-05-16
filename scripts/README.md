<!-- README.md is generated from README.Rmd. Please edit that file -->

# Scripts

## Observations

Data on the location of observed species of interest are taken from
published studies on myctophids, krill, salps and cephalopods and
combined using the [Process species presence
data.R](scripts/Process%20species%20presence%20data.R) script.

Random background points are generated using the `pseudoAbs.R`
[function](R/pseudoAbs.R) and combined using [Combine species presence
data with background
points.R](scripts/Combine%20species%20presence%20data%20with%20background%20points.R)

## Environmental covariates

Contemporary environmental covariates were processed using [Process
contemporary
covariates.R](scripts/Process%20contemporary%20covariates.R) and then
matched to observation locations using [Append covariates to presence
absence
data.R](scripts/Append%20covariates%20to%20presence%20absence%20data.R).
This script also outputs a raster [stack](data/covariate_stack.rds) of
the contemporary covariates for use in spatial predictions.

## Spatial blocking

Data were then split into 10 cross-validation folds by dividing the
Southern Ocean into 500 km hexagons with the `spatialsample` R package
and randomly assigning each hexagon into 1 of 10 blocks with [Assign
folds with
spatialsample.R](scripts/Assign%20folds%20with%20spatialsample.R). The
output is a list with each element containing the 10
[folds](data/folds.rds) for each species.

## Covariate selection

The optimal number of environmental covariates to describe the habitat
of each species were calculated using [Species-specific covariate
selection.R](scripts/Species-specific%20covariate%20selection.R). The
output of this script is a list with each element containing the model
ranking table for each species.
