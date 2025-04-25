# point_density_simulation

This repository stores an example dataset and R codes to simulate different combinations of the number of photoquadrats and the number of annotation points per photoquadrat to estimate live coral cover through the benthic point annotation approach.\
\
The example data consists of annotation data from 9 different reef plots.\
\
The codes, for each iteration, generate different photoquadrats/annotation point combination. For example, if we are simulating 10, 20, 30 and 40 photoquadrats and 10, 20, 50, 100, 200, 400, 600, 800 and 1000 annotation points per photoquadrats, each iteration will generate 36 sets of simulate data (i.e. 1 under each combination). For each set of simulated data, the codes will perform 10,000 bootstrap resampling to calculate mean coral cover and 95% confidence interval around the mean under two scenarios: option 1 and option 2.\
\
To perform the simulation 1,000 times for each of the photoquadrats/annotation point combination, you will need to run the codes for a total of 1,000 iterations. You can edit the codes to run multiple iterations at a time, but we recommend to keep it small (e.g. up to 20) as each iteration takes some time. The default setting is to run 1 iteration at a time.
