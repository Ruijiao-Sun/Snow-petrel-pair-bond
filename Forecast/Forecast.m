clc;clear;
%%
%load climate covariates
ncdisp('CESM2/CESM2-LE-Dec-Jan-Feb-snowdays_1850-2100_v2.nc');
Dec_nsnow = ncread('CESM2/CESM2-LE-Dec-Jan-Feb-snowdays_1850-2100_v2.nc','Dec_snow_days');
Jan_nsnow = ncread('CESM2/CESM2-LE-Dec-Jan-Feb-snowdays_1850-2100_v2.nc','Jan_snow_days');
Feb_nsnow = ncread('CESM2/CESM2-LE-Dec-Jan-Feb-snowdays_1850-2100_v2.nc','Feb_snow_days');
nsnow_raw = Dec_nsnow(:,1:250)+Jan_nsnow(:,2:251)+Feb_nsnow(:,2:251); %season 1850-2099
ncdisp('CESM2/December_max_2m_air_temp_1850-2100_v2.nc');
Dec_t = ncread('CESM2/December_max_2m_air_temp_1850-2100_v2.nc','max_2m_Temp_Dec');
Jan_t = ncread('CESM2/January_max_2m_air_temp_1850-2100_v2.nc','max_2m_Temp_Jan');
Feb_t = ncread('CESM2/February_max_2m_air_temp_1850-2100_v2.nc','max_2m_Temp_Feb');
Dec_t = Dec_t';Jan_t=Jan_t';Feb_t=Feb_t';
maxT_raw = (Dec_t(:,1:250)+Jan_t(:,2:251)+Feb_t(:,2:251))./3; %season 1850-2099
mu_obs = mean(maxT_raw(:,115:163),2);
sd_obs = std(maxT_raw(:,115:163)');
maxT = (maxT_raw(:,115:end)-mu_obs)./sd_obs';

mu_obs = mean(nsnow_raw(:,115:163),2);
sd_obs = std(nsnow_raw(:,115:163)');
nsnow = (nsnow_raw(:,115:end)-mu_obs)./sd_obs';
 %standardize using z-score
clear Jan_nsnow; clear Dec_nsnow; clear Feb_nsnow; clear Dec_t; clear Feb_t; clear Jan_t;
ncdisp('CESM2_aice/data_aice_coupled_breed_snpt.nc');
aice_breed = ncread('CESM2_aice/data_aice_coupled_breed_snpt.nc','layer');
aice_breed=aice_breed';
mu_obs = mean(aice_breed(:,65:113),2);
sd_obs = std(aice_breed(:,65:113)');
aice_breed = (aice_breed(:,65:end)-mu_obs)./sd_obs';
aice_nonbreed = ncread('CESM2_aice/data_aice_coupled_nonbreed_snpt.nc','layer');
aice_nonbreed = aice_nonbreed';
mu_obs = mean(aice_nonbreed(:,65:113),2);
sd_obs = std(aice_nonbreed(:,65:113)');
aice_nonbreed = (aice_nonbreed(:,65:(end-1))-mu_obs)./sd_obs';
%% estimates mu-mean,sigma-standard deviation from JAGS model output
%survival, estimates obtained from SPaice_8_F.RData
mu_s= 0.948;
sigma_s = 0.0076;
mu_s_aice_nonbreed = 0.152;
sigma_s_aice_nonbreed= 0.140;

%widowhood, estimates obtained from SPaice_8_F.RData
mu_w=0.0472 ;
sigma_w=0.00487;
mu_w_aice=0.094;
sigma_w_aice=0.096;

%breeding success, estimates obtained from SP_bs_full_F.RData 
mu_bs = 0.596;
sigma_bs = 0.033;
mu_bs_nsnow =-0.133 ;
sigma_bs_nsnow = 0.120;
mu_bs_maxT = 0.187;
sigma_bs_maxT = 0.129;
mu_bs_aice = -0.019;
sigma_bs_aice = 0.136;
% 
%divorce, estimates obtained from SP_d_full_F.RData
 mu_d_fail = 0.0323;
sigma_d_fail = 0.005;
mu_d_success = 0.0164;
sigma_d_success = 0.003;
% 
mu_d_nsnow_s = 0.190;
sigma_d_nsnow_s = 0.129;
mu_d_maxT_s = 0.064;
sigma_d_maxT_s = 0.141;
mu_d_aice_s = 0.043;
sigma_d_aice_s = 0.155;
mu_d_nsnow_f = 0.001;
sigma_d_nsnow_f = 0.119;
mu_d_maxT_f = -0.083;
sigma_d_maxT_f = 0.141;
mu_d_aice_f = 0.046;
sigma_d_aice_f = 0.145;

%%
nsim=100;
nenv = 50; %50 climate model ensembles
nind = 20000; %number of simulated life-history for each ensemble
divorce_rate = zeros(nenv*nsim,136); %population divorce rate projection from 1965-2099

nsimu=1;
for ens=1:nenv %loop through each ensemble
    for sim=1:nsim % loop through demographic error
        %creat vector for the probability of success
        p_success = invlogit(logit(normrnd(mu_bs,sigma_bs,[1,136])) + normrnd(mu_bs_maxT,sigma_bs_maxT,[1,136]).*maxT(ens,:)+normrnd(mu_bs_nsnow,sigma_bs_nsnow,[1,136]).*nsnow(ens,:)+normrnd(mu_bs_aice,sigma_bs_aice,[1,136]).*aice_breed(ens,:));
        p_survival = invlogit(logit(normrnd(mu_s,sigma_s,[1,136])) + normrnd(mu_s_aice_nonbreed,sigma_s_aice_nonbreed,[1,136]).*aice_nonbreed(ens,:));
        int_divorce_s = normrnd(mu_d_success,sigma_d_success,[1,136]);
        b1divorce_s = normrnd(mu_d_nsnow_s,sigma_d_nsnow_s,[1,136]);
        b2divorce_s = normrnd(mu_d_maxT_s,sigma_d_maxT_s,[1,136]);
        b3divorce_s = normrnd(mu_d_aice_s,sigma_d_aice_s,[1,136]);
        int_divorce_f = normrnd(mu_d_fail,sigma_d_fail,[1,136]);
        b1divorce_f = normrnd(mu_d_nsnow_f,sigma_d_nsnow_f,[1,136]);
        b2divorce_f = normrnd(mu_d_maxT_f,sigma_d_maxT_f,[1,136]);
        b3divorce_f = normrnd(mu_d_aice_f,sigma_d_aice_f,[1,136]);
        %b1divorce=normrnd(mu_divorce_success,sigma_divorce_success,[1,136]);
        %b2divorce=normrnd(mu_divorce_nsnow,sigma_divorce_nsnow,[1,136]);
        p_widow=invlogit(logit(normrnd(mu_w,sigma_w,[1,136])) + normrnd(mu_w_aice,sigma_w_aice,[1,136]).*aice_nonbreed(ens,:));

        sim_survival = zeros(nind,136);
         sim_success = zeros(nind,136); %creat matrice to store simulated breeding success
         sim_divorce = zeros(nind,136); %creat matrice to store simulated divorce
           sim_widow = zeros(nind,136); 

        for i=1:nind %loop through each life-history
            sim_survival(i,:)=binornd(1,p_survival);
            sim_success(i,:)=binornd(1,p_success); %simulated success
            sim_widow(i,:) = binornd(1,p_widow);
            %create vector for the probability of divorce depending on success
            p_divorce= invlogit(logit(int_divorce_s) + b1divorce_s.*nsnow(ens,:)+ b2divorce_s.*maxT(ens,:) + b3divorce_s.*aice_breed(ens,:)).*sim_success(i,:) + ...
                       invlogit(logit(int_divorce_f) + b1divorce_f.*nsnow(ens,:)+ b2divorce_f.*maxT(ens,:) + b3divorce_f.*aice_breed(ens,:)).*(1-sim_success(i,:));
            sim_divorce(i,:)=binornd(1,p_divorce); %simulated divorce
        end
        survival_rate(nsimu,:)=sum(sim_survival,1)./nind;
        divorce_rate(nsimu,:)=sum(sim_divorce,1)./nind; 
        success_rate(nsimu,:)=sum(sim_success,1)./nind;
        widow_rate(nsimu,:)=sum(sim_widow,1)./nind;
        nsimu=nsimu+1;
    end
end
%% store simulation for plotting
save("survival_rate_female.mat","survival_rate");
save("success_rate_female.mat","success_rate");
save("widow_rate_female.mat","widow_rate");
save("divorce_rate_female.mat","divorce_rate");
%%
%male
%survival, estimate obtained from SPaice_8_M.RData
mu_s= 0.948;
sigma_s = 0.007;
mu_s_aice_nonbreed = -0.0648;
sigma_s_aice_nonbreed= 0.139;

%widowhood, estimate obtained from SPaice_8_M.RData
mu_w=0.049;
sigma_w=0.007;
mu_w_aice=-0.040;
sigma_w_aice=0.130;

%breeding success, estimate obtained from SP_bs_full_M.RData
mu_bs = 0.589;
sigma_bs =0.037 ;
mu_bs_nsnow = -0.122;
sigma_bs_nsnow = 0.148;
mu_bs_maxT = 0.259;
sigma_bs_maxT =0.150;
mu_bs_aice = 0.052;
sigma_bs_aice =0.153 ;

%divorce, estimates obtained from SP_d_full_M.RData
mu_d_fail = 0.0441; %intercept
sigma_d_fail =0.0058 ;
mu_d_success =0.0184;
sigma_d_success = 0.003;

mu_d_nsnow_s =0.294;
sigma_d_nsnow_s =0.124 ;
mu_d_maxT_s =-0.033 ;
sigma_d_maxT_s =0.137 ;
mu_d_aice_s = -0.003;
sigma_d_aice_s =0.153 ;
mu_d_nsnow_f = 0.083;
sigma_d_nsnow_f = 0.104;
mu_d_maxT_f = -0.095;
sigma_d_maxT_f =0.128 ;
mu_d_aice_f = -0.097;
sigma_d_aice_f = 0.130;



%%
nsim=100;
nenv = 50; %50 climate model ensembles
nind = 20000; %number of simulated life-history for each ensemble
divorce_rate = zeros(nenv*nsim,136); %population divorce rate projection from 1965-2099

nsimu=1;
for ens=1:nenv %loop through each ensemble
    for sim=1:nsim % loop through demographic error
        %creat vector for the probability of success
        p_success = invlogit(logit(normrnd(mu_bs,sigma_bs,[1,136])) + normrnd(mu_bs_maxT,sigma_bs_maxT,[1,136]).*maxT(ens,:)+normrnd(mu_bs_nsnow,sigma_bs_nsnow,[1,136]).*nsnow(ens,:)+normrnd(mu_bs_aice,sigma_bs_aice,[1,136]).*aice_breed(ens,:));
        p_survival = invlogit(logit(normrnd(mu_s,sigma_s,[1,136])) + normrnd(mu_s_aice_nonbreed,sigma_s_aice_nonbreed,[1,136]).*aice_nonbreed(ens,:));
        int_divorce_s = normrnd(mu_d_success,sigma_d_success,[1,136]);
        b1divorce_s = normrnd(mu_d_nsnow_s,sigma_d_nsnow_s,[1,136]);
        b2divorce_s = normrnd(mu_d_maxT_s,sigma_d_maxT_s,[1,136]);
        b3divorce_s = normrnd(mu_d_aice_s,sigma_d_aice_s,[1,136]);
        int_divorce_f = normrnd(mu_d_fail,sigma_d_fail,[1,136]);
        b1divorce_f = normrnd(mu_d_nsnow_f,sigma_d_nsnow_f,[1,136]);
        b2divorce_f = normrnd(mu_d_maxT_f,sigma_d_maxT_f,[1,136]);
        b3divorce_f = normrnd(mu_d_aice_f,sigma_d_aice_f,[1,136]);
        %b1divorce=normrnd(mu_divorce_success,sigma_divorce_success,[1,136]);
        %b2divorce=normrnd(mu_divorce_nsnow,sigma_divorce_nsnow,[1,136]);
        p_widow=invlogit(logit(normrnd(mu_w,sigma_w,[1,136])) + normrnd(mu_w_aice,sigma_w_aice,[1,136]).*aice_nonbreed(ens,:));

        sim_survival = zeros(nind,136);
         sim_success = zeros(nind,136); %creat matrice to store simulated breeding success
         sim_divorce = zeros(nind,136); %creat matrice to store simulated divorce
           sim_widow = zeros(nind,136); 

        for i=1:nind %loop through each life-history
            sim_survival(i,:)=binornd(1,p_survival);
            sim_success(i,:)=binornd(1,p_success); %simulated success
            sim_widow(i,:) = binornd(1,p_widow);
            %create vector for the probability of divorce depending on success
            p_divorce= invlogit(logit(int_divorce_s) + b1divorce_s.*nsnow(ens,:)+ b2divorce_s.*maxT(ens,:) + b3divorce_s.*aice_breed(ens,:)).*sim_success(i,:) + ...
                       invlogit(logit(int_divorce_f) + b1divorce_f.*nsnow(ens,:)+ b2divorce_f.*maxT(ens,:) + b3divorce_f.*aice_breed(ens,:)).*(1-sim_success(i,:));
            sim_divorce(i,:)=binornd(1,p_divorce); %simulated divorce
        end
        survival_rate_male(nsimu,:)=sum(sim_survival,1)./nind;
        divorce_rate_male(nsimu,:)=sum(sim_divorce,1)./nind; 
        success_rate_male(nsimu,:)=sum(sim_success,1)./nind;
        widow_rate_male(nsimu,:)=sum(sim_widow,1)./nind;
        nsimu=nsimu+1;
    end
end
%% store simulation for plotting
save("survival_rate_male.mat","survival_rate_male");
save("success_rate_male.mat","success_rate_male");
save("widow_rate_male.mat","widow_rate_male");
save("divorce_rate_male.mat","divorce_rate_male");
