load("SP0bs_F_lower.mat")
load("SP0bs_F_mean.mat")
load("SP0bs_F_upper.mat")
load("SP0bs_M_lower.mat")
load("SP0bs_M_mean.mat")
load("SP0bs_M_upper.mat")
fig=figure;
fig.Position=[560,240,350,400]
%xtick location
f_x = 0.9; hold on
failed_f_x = f_x+0.6;
m_x = f_x+0.1; 
failed_m_x = m_x+0.6;
%estimate
f_y = F_mean.d_SB;
f_yneg = F_lower.d_SB-f_y;
f_ypos = F_upper.d_SB-f_y;
%individuals failed at their previous breeding attempts
failed_f_y = F_mean.d_SBf;
failed_f_yneg = F_lower.d_SBf-failed_f_y;
failed_f_ypos = F_upper.d_SBf-failed_f_y;
m_y = M_mean.d_SB;
m_yneg = M_lower.d_SB-m_y;
m_ypos = M_upper.d_SB-m_y;
%individuals failed at their previous breeding attempts
failed_m_y = M_mean.d_SBf;
failed_m_yneg = M_lower.d_SBf-failed_m_y;
failed_m_ypos = M_upper.d_SBf-failed_m_y;

errorbar(f_x,f_y,f_yneg,f_ypos,'o','Color',[153 51 51]/256,'MarkerEdgeColor','k','MarkerFaceColor',[239 37 119]/256,"Markersize",10);
errorbar(failed_f_x,failed_f_y,failed_f_yneg,failed_f_ypos,'o','Color',[153 51 51]/256,'MarkerEdgeColor','k','MarkerFaceColor',[239 37 119]/256,"Markersize",10);
errorbar(m_x,m_y,m_yneg,m_ypos,'o','Color',[0 102 204]/256,'MarkerEdgeColor','k','MarkerFaceColor',[11 141 221]/256,"Markersize",10);
errorbar(failed_m_x,failed_m_y,failed_m_yneg,failed_m_ypos,'o','Color',[0 102 204]/256,'MarkerEdgeColor','k','MarkerFaceColor',[11 141 221]/256,"Markersize",10);
xlim([0.5 2])
ylim([0 0.1])
xticks([f_x+0.05,f_x+0.6])
xticklabels({'Success', 'Failure', 'ND','NB'})
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica')
pbaspect([1.2 1.2 1]);
box on
%legend("F","M",'Location',"southeast");
set(gca,'linewidth',1.1,'fontsize',14,'fontname','Helvetica')
ylabel('Divorce rate');xlabel('Previous breeding outcomes');
%print(gcf,'divorce.pdf','-dpdf','-r600'); 
exportgraphics(fig,'divorce_mean.pdf','BackgroundColor','none')
exportgraphics(fig,'divorce_mean.png','BackgroundColor','none')