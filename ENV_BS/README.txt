File list: 

SP_pair.RData contains individual life histories as a series of life events, including survival, widowhood, divorce, breeding, and breeding success (details provided below). These life histories are then linked to detection data. In our study, all individuals enter the model at their second breeding attempt during adulthood, as we needed to determine whether they were breeding with the same partner or a new one. The raw observation data spans from 1963 to 2016, while the life history data for estimates ends in 2012 to determine the survival status for individuals still alive in 2012.
SURVIVAL: Indicates whether an individual was alive or not. "1" – alive, "NA" – no information available.
WIDOWHOOD: Indicates whether an individual was widowed or not. "1" – widowed, "0" – not widowed.
DIVORCE: Indicates whether an individual was divorced or not. "1" – divorced, "0" – not divorced.
BREEDING: Indicates whether an individual bred or not. "1" – bred, "0" – did not breed, "NA" – no information available.
SUCCESS: Indicates whether an individual bred successfully or not, conditional on breeding. "1" – succeeded, "0" – failed, "NA" – no information available.
DETECTED: Indicates whether an individual was detected or not. "1" – detected, "0" – not detected.
BS: Breeding success of the previous breeding attempt. "1" – previously successful, "0" – previously failed, "NA" – no information available.
SEX: Sex of an individual. "F" – Female, "M" – Male.
year_B: Year of the previous breeding attempt.

SP_3_full.txt contains a JAGS model to test the effects of environmental conditions during breeding season on breeding success.

JAGS_SP_3F_full.R runs jags model SP_3_full.txt for females through the R package jagsUI using data SP_pair.RData.
JAGS_SP_3M_full.R runs jags model SP_3_full.txt for males through the R package jagsUI using data SP_pair.RData.

Fig4efg.m is a Matlab code to plot Fig. 4efg in the manuscript main text using model output. 

aice_breed.1964.2012.xlsx contains sea aice area covariate.
mean_maxT.breed.1964.2012.xlsx contains mean daily maximum temperature covariate.
n_snow.breed.1964.2012.xlsx contains number of snow days covariate.

dataout.R extract the posterior from model output.

Model output files include: 
posterior_SP_bs_ENV.mat
SP_bs_full_F.RData
SP_bs_full_M.RData
