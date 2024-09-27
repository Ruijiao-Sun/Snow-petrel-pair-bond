%% 
clc;clear;
nsim=100;
nenv = 50; %50 climate model ensembles
load("divorce_rate_female.mat");
load("divorce_rate_male.mat");
load("widow_rate_male.mat");
load("widow_rate_female.mat");
load("survival_rate_female.mat");
load("survival_rate_male.mat");
load("success_rate_male.mat");
load("success_rate_female.mat");
%%
rand_i = randi(5000,[5 1]);
mean_d = mean(divorce_rate_male,1);
sorted_d=sort(divorce_rate_male);
lower=sorted_d(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_d(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[0 102 204]/256)
plot(pro_year,upper_pro,"--","Color",[0 102 204]/256)
ciplot( lower_pro,upper_pro,pro_year ,[153, 204, 255]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,divorce_rate_male(rand_i,49:end))
plot(2012:2099,mean_d(49:end),"LineWidth",2,"Color",[0 102 204]/256);
ylim([0 0.25]);
xlim([2010 2100]);
xlabel("Year");ylabel("Male divorce rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out')
exportgraphics(fig,'figure\M_divorce_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\M_divorce_forecast.png','BackgroundColor','none')
%%
rand_i = randi(5000,[5 1]);
mean_w = mean(widow_rate_male,1);
sorted_w=sort(widow_rate_male);
lower=sorted_w(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_w(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[0 102 204]/256)
plot(pro_year,upper_pro,"--","Color",[0 102 204]/256)
ciplot( lower_pro,upper_pro,pro_year ,[153, 204, 255]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,widow_rate_male(rand_i,49:end))
plot(2012:2099,mean_w(49:end),"LineWidth",2,"Color",[0 102 204]/256);
ylim([0 0.25]);
xlim([2010 2100]);
xlabel("Year");ylabel("Male widowhood rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off");
exportgraphics(fig,'figure\M_widow_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\M_widow_forecast.png','BackgroundColor','none')
%%
rand_i = randi(5000,[5 1]);
mean_d = mean(divorce_rate,1);
sorted_d=sort(divorce_rate);
lower=sorted_d(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_d(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[153 51 51]/256)
plot(pro_year,upper_pro,"--","Color",[153 51 51]/256)
ciplot( lower_pro,upper_pro,pro_year,[225 175 175]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,divorce_rate(rand_i,49:end))
plot(2012:2099,mean_d(49:end),"LineWidth",2,"Color",[153 51 51]/256);
ylim([0 0.25]);
xlim([2010 2100]);
xlabel("Year");ylabel("Female divorce rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out')
exportgraphics(fig,'figure\F_divorce_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\F_divorce_forecast.png','BackgroundColor','none')
%%
rand_i = randi(5000,[5 1]);
mean_w = mean(widow_rate,1);
sorted_w=sort(widow_rate);
lower=sorted_w(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_w(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[153 51 51]/256)
plot(pro_year,upper_pro,"--","Color",[153 51 51]/256)
ciplot( lower_pro,upper_pro,pro_year ,[225 175 175]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,widow_rate(rand_i,49:end))
plot(2012:2099,mean_w(49:end),"LineWidth",2,"Color",[153 51 51]/256);
ylim([0 0.25]);
xlim([2010 2100]);
xlabel("Year");ylabel("Female widowhood rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off")
exportgraphics(fig,'figure\F_widow_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\F_widow_forecast.png','BackgroundColor','none')

%% F survival
rand_i = randi(5000,[5 1]);
mean_s = mean(survival_rate,1);
sorted_s=sort(survival_rate);
lower=sorted_s(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_s(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[153 51 51]/256)
plot(pro_year,upper_pro,"--","Color",[153 51 51]/256)
ciplot( lower_pro,upper_pro,pro_year ,[225 175 175]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,survival_rate(rand_i,49:end))
plot(2012:2099,mean_s(49:end),"LineWidth",2,"Color",[153 51 51]/256);
ylim([0.4 1]);
xlim([2010 2100]);
xlabel("Year");ylabel("Female survival rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off")
exportgraphics(fig,'figure\F_survival_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\F_survival_forecast.png','BackgroundColor','none')
%%
rand_i = randi(5000,[5 1]);
mean_s = mean(survival_rate_male,1);
sorted_s=sort(survival_rate_male);
lower=sorted_s(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_s(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[0 102 204]/256)
plot(pro_year,upper_pro,"--","Color",[0 102 204]/256)
ciplot( lower_pro,upper_pro,pro_year ,[153, 204, 255]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,survival_rate_male(rand_i,49:end))
plot(2012:2099,mean_s(49:end),"LineWidth",2,"Color",[0 102 204]/256);
ylim([0.4 1]);
xlim([2010 2100]);
xlabel("Year");ylabel("Male survival rate");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off");
exportgraphics(fig,'figure\M_survival_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\M_survival_forecast.png','BackgroundColor','none')
%%
%% F success
rand_i = randi(5000,[5 1]);
mean_bs = mean(success_rate,1);
sorted_bs=sort(success_rate);
lower=sorted_bs(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_bs(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[153 51 51]/256)
plot(pro_year,upper_pro,"--","Color",[153 51 51]/256)
ciplot( lower_pro,upper_pro,pro_year ,[225 175 175]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,success_rate(rand_i,49:end))
plot(2012:2099,mean_bs(49:end),"LineWidth",2,"Color",[153 51 51]/256);
ylim([0 1]);
xlim([2010 2100]);
xlabel("Year");ylabel("Female breeding success");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off")
exportgraphics(fig,'figure\F_success_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\F_success_forecast.png','BackgroundColor','none')
%% male success
rand_i = randi(5000,[5 1]);
mean_bs = mean(success_rate_male,1);
sorted_bs=sort(success_rate_male);
lower=sorted_bs(nenv*nsim*0.025,:); lower_obs = lower(1:49); lower_pro = lower(49:end);
upper=sorted_bs(nenv*nsim*0.975,:);  upper_obs = upper(1:49); upper_pro = upper(49:end);
obs_year=1964:2012; %year of observation
pro_year=2012:2099; %year of project
fig=figure;
fig.Position=[560,240,350,400]
hold on;
pbaspect([1 1 1]);
%plot(1964:2099,divorce_rate,"Color",[0,0,0,0.01]);
%plot(1964:2012,lower_obs,"--","Color","k")
%plot(1964:2012,upper_obs,"--","Color","k")
%ciplot( lower_obs,upper_obs,obs_year ,[0.75, 0.8, 0.9]);
%plot(obs_year,mean_d(1:49),"LineWidth",1.5,"Color","k");
plot(pro_year,lower_pro,"--","Color",[0 102 204]/256)
plot(pro_year,upper_pro,"--","Color",[0 102 204]/256)
ciplot(lower_pro,upper_pro,pro_year ,[153, 204, 255]/256);
colororder([255/256 204/256 102/256; 255/256 153/256 204/256; 0 153/256 153/256;179/256 179/256 255/256;140/256 140/256 140/256]);
plot(2012:2099,success_rate_male(rand_i,49:end))
plot(2012:2099,mean_bs(49:end),"LineWidth",2,"Color",[0 102 204]/256);
ylim([0 1]);
xlim([2010 2100]);
xlabel("Year");ylabel("Male breeding success");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica','TickDir','out',"Clipping","off");
exportgraphics(fig,'figure\M_success_forecast.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\M_success_forecast.png','BackgroundColor','none')

%%
trend_fw=[];
for i = 1:5000
    p = polyfit(1:136,widow_rate(i,:),1);
    p=p(1);
    trend_fw = [trend_fw,p];
end
sum(trend_fw>0)/5000
fig=figure
fig.Position=[560,240,150,150]
h=histogram(trend_fw,25,'Normalization','probability','EdgeAlpha',0.5,'FaceAlpha',0.4,"FaceColor",[0.8500 0.3250 0.0980]);hold on;box off; xlim([-0.3e-3 0.3e-3])
xlabel('Trend'); ylabel('Density')
ax = gca;                        % gets the current axes
ax.XAxisLocation = 'origin';     % sets them to zero
ax.YAxisLocation = 'origin';     % sets them to zero
ax.Box = 'off';             
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica','TickDir','out');
exportgraphics(fig,'figure\Fw_trend.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\Fw_trend_forecast.png','BackgroundColor','none')
%100% negative
%%
trend_mw=[];
for i = 1:5000
    p = polyfit(1:136,widow_rate_male(i,:),1);
    p=p(1);
    trend_mw = [trend_mw,p];
end
sum(trend_mw>0)/5000
fig=figure
fig.Position=[560,240,400,400]
plot(2012:2099,widow_rate_male(120,end-87:end),LineWidth=1.5);hold on;box off;xlim([2010 2100])
p=polyfit(1:136,widow_rate_male(120,:),1);
plot(2012:2099,0.0489+0.0002.*((136-87):136),LineWidth=1.5);
legend("Trajectory","Regression trend");ylabel("Widowhood probability");xlabel("Year")
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica','TickDir','out');
exportgraphics(fig,'figure\Mw_trend_example.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\Mw_trend_example.png','BackgroundColor','none')
%%
fig=figure
fig.Position=[560,240,150,150]
h=histogram(trend_mw,'Normalization','probability','EdgeAlpha',0.5,'FaceAlpha',0.4,"FaceColor",[0 0.4470 0.7410]);hold on;box off; xlim([-0.7e-3 0.7e-3])
xlabel('Trend'); ylabel('Density')
ax = gca;                        % gets the current axes
ax.XAxisLocation = 'origin';     % sets them to zero
ax.YAxisLocation = 'origin';     % sets them to zero
ax.Box = 'off';    
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica','TickDir','out');
exportgraphics(fig,'figure\Mw_trend.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\Mw_trend_forecast.png','BackgroundColor','none')
%100% positive
%58%
%%
trend_fd=[];
for i = 1:5000
    p = polyfit(1:136,divorce_rate(i,:),1);
    p=p(1);
    trend_fd = [trend_fd,p];
end
%sum(trend_fd>0)/50
fig=figure
fig.Position=[560,240,150,150]
h=histogram(trend_fd,'Normalization','probability','EdgeAlpha',0.5,'FaceAlpha',0.4,"FaceColor",[0.8500 0.3250 0.0980]);hold on;box off; xlim([-0.25e-3 0.25e-3])
xlabel('Trend'); ylabel('Density')
ax = gca;                        % gets the current axes
ax.XAxisLocation = 'origin';     % sets them to zero
ax.YAxisLocation = 'origin';     % sets them to zero
ax.Box = 'off';             
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica','TickDir','out');
exportgraphics(fig,'figure\Fd_trend.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\Fd_trend_forecast.png','BackgroundColor','none')

%46% 
%%
trend_md=[];
for i = 1:5000
    p = polyfit(1:136,divorce_rate_male(i,:),1);
    p=p(1);
    trend_md = [trend_md,p];
end
sum(trend_md>0)/5000
fig=figure
fig.Position=[560,240,150,150]
h=histogram(trend_md,'Normalization','probability','EdgeAlpha',0.5,'FaceAlpha',0.4,"FaceColor",[0 0.4470 0.7410]);hold on;box off; xlim([-0.25e-3 0.25e-3])
xlabel('Trend'); ylabel('Density')
ax = gca;                        % gets the current axes
ax.XAxisLocation = 'origin';     % sets them to zero
ax.YAxisLocation = 'origin';     % sets them to zero
ax.Box = 'off';    
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica','TickDir','out');
exportgraphics(fig,'figure\Md_trend.pdf','BackgroundColor','none')
exportgraphics(fig,'figure\Md_trend_forecast.png','BackgroundColor','none')

% 54%
%%
%%
trend_fs=[];
for i = 1:50
    p = polyfit(1:136,survival_rate(i,:),1);
    p=p(1);
    trend_fs = [trend_fs,p];
end
sum(trend_fs>0)/50
%0% 

