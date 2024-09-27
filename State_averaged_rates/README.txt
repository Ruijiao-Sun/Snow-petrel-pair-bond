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

SP_0.txt contains baseline JAGS model to estimate vital rates and pair-bond disruption probability that are not dependent on individual states. 

JAGS_SP_0F.R runs jags model SP_0.txt for females through the R package jagsUI using data SP_pair.RData.
JAGS_SP_0M.R runs jags model SP_0.txt for males through the R package jagsUI using data SP_pair.RData.

Fig3b.m is a Matlab code to plot Fig. 3b in the manuscript main text using model output. 

Model output files include: 
SP0bs_F_lower.matSP0bs_F_mean.matSP0bs_F_upper.matSP0bs_M_lower.matSP0bs_M_mean.matSP0bs_M_upper.matSPbs_0_F.RDataSPbs_0_M.RData