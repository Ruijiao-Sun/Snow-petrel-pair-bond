rm(list = ls())#clear lists
my.packs <- c("coda", "mvtnorm", "rjags", "dplyr", "reshape", "readxl",'jagsUI',"R.matlab")
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE,configure.args="--enable-rpath")
lapply(my.packs, require, character.only = TRUE)

setwd(getwd())
load("SP_pair.RData")
#cov1: nsnow
#cov2: meanmaxT
#cov3: aice
nsnow <- readxl::read_xlsx("n_snow.breed.1964.2012.xlsx",col_names = F)
nsnow <-as.matrix(nsnow)
COV1 <-scale(nsnow)

maxT <- readxl::read_xlsx("mean_maxT.breed.1964.2012.xlsx",col_names = F)
maxT <-as.matrix(maxT)
COV2 <-scale(maxT)

aice <- readxl::read_xlsx("aice_breed.1964.2012.xlsx",col_names = F)
aice <-as.matrix(aice)
COV3 <-scale(aice)

#range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))} #define a function to scale data between 0 and 1
#SAM<-range01(SAM,na.rm=T)
cov.cov1 = seq(min(COV1,na.rm=T),max(COV1,na.rm=T), length.out = 1000)
n.pred = length(cov.cov1)

cov.cov2 = seq(min(COV2,na.rm=T),max(COV2,na.rm=T), length.out = 1000)
n.pred = length(cov.cov2)

cov.cov3 = seq(min(COV3,na.rm=T),max(COV3,na.rm=T), length.out = 1000)
n.pred = length(cov.cov3)


#only select female individuals
DETECTED<-DETECTED[which(SEX=="F"),]
SURVIVAL<-SURVIVAL[which(SEX=="F"),]
BREEDING<-BREEDING[which(SEX=="F"),]
WIDOWHOOD<-WIDOWHOOD[which(SEX=="F"),]
DIVORCE <- DIVORCE[which(SEX=="F"),]
BS <- BS[which(SEX=="F"),]
SUCCESS <- SUCCESS[which(SEX=="F"),]
S <- S[which(SEX=="F"),]
year_B <- year_B[which(SEX=="F"),]

### for jags
### bundle data
cmr_data <- list(DETECTED = DETECTED,
                 SURVIVAL = SURVIVAL,
                 BREEDING = BREEDING,
                 DIVORCE = DIVORCE,
                 WIDOWHOOD = WIDOWHOOD,
                 SUCCESS = SUCCESS,
                 BS = BS,
                 S = S,
                 year_B = year_B); 
rm(SURVIVAL, BREEDING, DIVORCE, WIDOWHOOD, DETECTED,S,SUCCESS,BS)
jagsdata <- vector(mode = "list", length = 0)
jagsdata$n_year <- ncol(cmr_data$DETECTED)
jagsdata$n_ind <- nrow(cmr_data$DETECTED)
jagsdata$FIRST <- apply(cmr_data$DETECTED, 1, function(row_i) { min(which(row_i == 1)) })
jagsdata$LAST <- apply(cmr_data$SURVIVAL, 1, function(row_i) { ifelse(any(row_i == 0, na.rm = TRUE), min(which(row_i == 0)), ncol(cmr_data$SURVIVAL)) })
jagsdata$SURVIVAL <- cmr_data$SURVIVAL
jagsdata$BREEDING <- cmr_data$BREEDING
jagsdata$WIDOWHOOD <- cmr_data$WIDOWHOOD
jagsdata$DIVORCE <- cmr_data$DIVORCE
jagsdata$DETECTED <- cmr_data$DETECTED
jagsdata$S <- cmr_data$S
jagsdata$BS <- cmr_data$BS
jagsdata$SUCCESS <- cmr_data$SUCCESS
jagsdata$year_B <- cmr_data$year_B
### function to initialize sampler

initfct <- function(idchain, jagsdata, perfect_detection = TRUE) 
{
  ### covariances
  tauA <- rgamma(4, shape = 1.5, rate = 1.5)
  A <- LA <- matrix(NA, nrow = 4, ncol = 4)
  diag(A) <- abs(rnorm(4)) * 1.5
  LA[2, 1] <- 0.25 * rnorm(1)
  LA[3, 1] <- 0.25 * rnorm(1)
  LA[4, 1] <- 0.25 * rnorm(1)
  LA[3, 2] <- 0.25 * rnorm(1)
  LA[4, 2] <- 0.25 * rnorm(1)
  LA[4, 3] <- 0.25 * rnorm(1)

  tauE <- rgamma(4, shape = 1.5, rate = 1.5)
  E <- LE <- matrix(NA, nrow = 4, ncol = 4)
  diag(E) <- abs(rnorm(4)) * 1.5
  LE[2, 1] <- 0.25 * rnorm(1)
  LE[3, 1] <- 0.25 * rnorm(1)
  LE[4, 1] <- 0.25 * rnorm(1)
  LE[3, 2] <- 0.25 * rnorm(1)
  LE[4, 2] <- 0.25 * rnorm(1)
  LE[4, 3] <- 0.25 * rnorm(1)

  init <- list(mu_phi = 1.5 * rnorm(1),
               mu_rho = 1.5 * rnorm(1),
               mu_psi = 1.5 * rnorm(1),
               mu_pi = 1.5 * rnorm(1),
               gamma = rnorm(4) * log(2) / 4,
               tauA = tauA,
               A = A,
               LA = LA,
               xi_a = replicate(4, rnorm(jagsdata$n_ind)) %*% diag(1 / sqrt(tauA)),
               tauE = tauE,
               E = E,
               LE = LE,
               xi_e = replicate(4, rnorm(jagsdata$n_year)) %*% diag(1 / sqrt(tauE))
               )
  
  if(!perfect_detection) 
  {
    ### generate coherent values for 1 when not detected
    widowhood <- divorce <- breed <- survival<- success <- with(jagsdata, array(NA, dim = c(n_ind, n_year)))
    # loop over individuals
    for(i in 1:jagsdata$n_ind) 
    {
      first_det <- jagsdata$FIRST[i]
      last_det <- max(which(jagsdata$SURVIVAL[i, ] == 1))
      #########survival first###########
      if(last_det + 1 < jagsdata$n_year) 
      {
        death <- sample((last_det + 1):jagsdata$n_year, size = 1)
      }
      else 
      {
        death <- jagsdata$n_year
      }
      survival[i, first_det:death] <- 1
      if(death != jagsdata$n_year) 
      {
        survival[i, (death + 1):jagsdata$n_year] <- 0
      }
      #########breeding#################
      if(last_det < jagsdata$n_year) 
      { # check that you need to do something
        if(any(jagsdata$DETECTED[i, ] == 0, na.rm = TRUE)) 
        {          
            #breed[i, first_det] <- 1
            #breed[i, (first_det+1):jagsdata$n_year] <- 0#rbinom(length((last_det+1):jagsdata$n_year), size = 1, prob = 0.5)
            bond_disruption_year <- sample((last_det + 1):jagsdata$n_year, size = 1)
            toss<-rbinom(1,size=1,prob=0.5)
            if (toss==1)
            {
             widowhood[i,(last_det + 1):jagsdata$n_year]<-0
             widowhood[i,bond_disruption_year]<-1
             divorce[i,(last_det + 1):jagsdata$n_year]<-0
            }
            else
            {
             widowhood[i,(last_det + 1):jagsdata$n_year]<-0
             divorce[i,(last_det + 1):jagsdata$n_year]<-0
             divorce[i,bond_disruption_year]<-1
            }
        }
      }
    }
  
     survival <- survival 
        breed <- breed * survival
        success <- success * breed
     widowhood <- widowhood * (1-divorce) * survival
       divorce <- divorce * (1-widowhood) * survival

    init$SURVIVAL <- survival * ifelse(is.na(jagsdata$SURVIVAL), 1, NA)
     init$BREEDING <- breed * ifelse(is.na(jagsdata$BREEDING), 1, NA)
     init$SUCCESS <- success * ifelse(is.na(jagsdata$SUCCESS), 1, NA)
    init$WIDOWHOOD <- widowhood * ifelse(is.na(jagsdata$WIDOWHOOD), 1, NA)
      init$DIVORCE <- divorce * ifelse(is.na(jagsdata$DIVORCE), 1, NA)
    init$mu_p <- 1.5 * rnorm(2)
    init$unscaled_sigma2_p <- rgamma(2, 1.0, 1.0)
    init$rho_p <- rgamma(2, 1.0, 1.0)
    
    # Not fixing detection
    init$upsilon <- cbind(init$mu_p[1] + 0.05 * cumsum(rnorm(jagsdata$n_year)),
                           init$mu_p[2] + 0.05 * cumsum(rnorm(jagsdata$n_year)),
                          init$mu_p[3] + 0.05 * cumsum(rnorm(jagsdata$n_year))
    )
   }
  return(init)
  }


#set initial condition
IV <- initfct(jagsdata = jagsdata, perfect_detection = F)  

library(jagsUI)  # The R package which makes the interface with JAGS
 
# Bundle data
jags.data <- list( SURVIVAL = jagsdata$SURVIVAL, 
                   BREEDING = jagsdata$BREEDING, 
                  WIDOWHOOD = jagsdata$WIDOWHOOD, 
                    DIVORCE = jagsdata$DIVORCE, 
                   DETECTED = jagsdata$DETECTED,
                          S = jagsdata$S,
                     year_B = jagsdata$year_B,
                    SUCCESS = jagsdata$SUCCESS,
                         BS = jagsdata$BS,
                      FIRST = jagsdata$FIRST, 
                       LAST = jagsdata$LAST, 
                      n_ind = jagsdata$n_ind, 
                     n_year = jagsdata$n_year,
                       COV1 = COV1,
                       COV2 = COV2,
                       COV3 = COV3,
                     n.pred = n.pred,
                    cov.cov1 = cov.cov1,
                    cov.cov2 = cov.cov2,
                    cov.cov3 = cov.cov3)

inits <- function(){list(SURVIVAL = IV$SURVIVAL, 
                         BREEDING = IV$BREEDING, 
                         WIDOWHOOD = IV$WIDOWHOOD,
                         DIVORCE = IV$DIVORCE)}

# Parameters monitored
parameters <- c(#female mean vital rates
                "s_SB", "s_NW", "s_ND", "s_NB", "s_D", "s_W",
                "s_SB_mean", "s_NW_mean", "s_ND_mean", "s_NB_mean", "s_D_mean", "s_W_mean",
                "w_SB", "w_NW", "w_ND", "w_NB", "w_D", "w_W",
                "w_SB_mean", "w_NW_mean", "w_ND_mean", "w_NB_mean", "w_D_mean", "w_W_mean",
                "d_SB", "d_NW", "d_ND", "d_NB", "d_D", "d_W",
                "d_s_mean",
                "d_SBf", "d_NWf", "d_NDf", "d_NBf", "d_Df", "d_Wf",
                "d_f_mean", 
                "beta_SB", "beta_NW", "beta_ND", "beta_NB", "beta_D", "beta_W",
                "beta_SB_mean", "beta_NW_mean", "beta_ND_mean", "beta_NB_mean", "beta_D_mean", "beta_W_mean",
                "detect_SB", "detect_NW", "detect_ND", "detect_NB", "detect_D", "detect_W",
                "detect_SB_mean", "detect_NW_mean", "detect_ND_mean", "detect_NB_mean", "detect_D_mean", "detect_W_mean",
                "f_SB", "f_NW", "f_ND", "f_NB", "f_D", "f_W",
                "f_SB_mean", "f_NW_mean", "f_ND_mean", "f_NB_mean", "f_D_mean", "f_W_mean",
                "slope1s","slope2s","slope3s","slope1f","slope2f","slope3f"
                )
# MCMC settings
ni <- 16000      # number of iteration 
nt <- 2          # Thinning interval
nb <- 8000       # number of burnin
nc <- 3          # number of chains

# Call JAGS from R (BRT 1 min)
out1 <- jags(jags.data, inits, parameters, "SP_4_full.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

save(out1,file = "SP_d_full_F.RData")

writeMat("SP_d_full_F_mean.mat",F_mean=out1$mean)
writeMat("SP_d_full_F_lower.mat",F_lower=out1$q2.5)
writeMat("SP_d_full_F_upper.mat",F_upper=out1$q97.5)









