<!-- README.md is generated from README.Rmd. Please edit that file -->

# Scripts

## Observations

Data on the location of observed species of interest are taken from
published studies on myctophids, krill, salps and cephalopods and
combined using the [Process species presence
data.R](scripts/Process%20species%20presence%20data.R) script.

Random background points are generated using the `buffered_points.R`
[function](R/buffered_points.R) and combined using [combine presences
and absences.R](scripts/combine%20presences%20and%20absences.R)

Contemporary environmental covariates were processed using [Process
contemporary covariates.R](R/Process%20contemporary%20covariates.R) and
then matched to observation locations using [Append covariates to
presence absence
data.R](Append%20covariates%20to%20presence%20absence%20data.R)
