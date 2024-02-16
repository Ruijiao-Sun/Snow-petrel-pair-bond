# Snow-petrel-pair-bond
Temporary data repository for manuscript submitted to Ecology Letters during initial submission entitled "The Dynamics of Pair Bonds and Their Demographic
Consequences in Variable Environments and a Shifting Climate"

2024-02-16

Software:
These files includes pair-bond observation events to estimate vital rates and probabilities of divorce and widowhood using JAGS through R package JagsUI. 

Author(s)
Ruijiao Sun
Woods Hole Oceanographic Institution
Massachusetts Institute of Technology
University of California, Santa Barbara
rsun@whoi.edu
ruijiaos@ucsb.edu

File list:
SP_pair.txt
SP_JAGS.txt

Description:

These data SP_pair.txt come from long-term monitoring of snow petrels breed at Ile des Pétrels, Pointe Géologie Archipelago (66°40' S, 140°01' E), Terre Adélie, Antarctica. Since 1963, an annual long-term monitoring study has been conducted. Adults and chicks were leg-banded with stainless-steel bands, with adult sex determined by vocalization and relative size. Nest surveys during incubation and fledging periods determined breeding success and pair identities.

The pair-bond state code corresponded to the states defined in the main text Fig 2. "0" indicates no detection and thus no information in the current year. And lowercase indicates no detection. In our study, all individuals entered the model at their second breeding attempt in the adult stage as we needed to define whether they were breeding with the same previous partner or not. 

SP_JAGS.txt included the JAGS code to run the baseline model.
