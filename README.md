# categorical_attraction_SCR
R script for running MCMC for a spatial capture recapture model with attractions between individuals.

This R script fits the spatial capture-recapture model with categorical attractions between individuals from the paper:

McLaughlin, P., & Bar, H. (2021). A spatial capture–recapture model with attractions between individuals. Environmetrics, 32(1), e2653.

The model is fit via MCMC using the Metropolis–Hastings algorithm. Users must specify a capture history data.frame at the start of the script which consists of individual IDs, where and when each individual was captured, as well as coordinates for trapping locations. Sample input data is given at the start of the script.
