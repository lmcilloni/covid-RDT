clear all;
% load('dynamic_results_ME_220418.mat');
% load('dynamic_results_ME_RDT_220118.mat');
load('dynamic_results_ME_221121.mat')
full_inc_ME = full_inc;
pcr_cons_ME = pcr_cons;
lfa_cons_ME = lfa_cons;
fp_quar_ME  = fp_quar;
ca_ME       = ca;
% load('dynamic_results_LO_220418.mat');
% load('dynamic_results_LO_RDT_220118.mat');
load('dynamic_results_LO_221121.mat')
full_inc_LO = full_inc;
pcr_cons_LO = pcr_cons;
lfa_cons_LO = lfa_cons;
fp_quar_LO  = fp_quar;
ca_LO       = ca;
% load('dynamic_results_HI_220418.mat');
% load('dynamic_results_HI_RDT_220118.mat');
load('dynamic_results_HI_221121.mat')
full_inc_HI = full_inc;
pcr_cons_HI = pcr_cons;
lfa_cons_HI = lfa_cons;
fp_quar_HI  = fp_quar;
ca_HI       = ca;


figure;
pl1=plot(1:tf, mean(full_inc_ME),'linewidth',2.5);hold on;
% errorbar(find(mean(full_inc_ME)==max(mean(full_inc_ME))), max(mean(full_inc_ME)), ...
%     max(mean(full_inc_ME))-max(quantile(full_inc_ME,0.025)), max(quantile(full_inc_ME,0.975))-max(mean(full_inc_ME)), 'color',[0.00,0.45,0.74])
%jbfill(1:tf,quantile(full_inc_ME,0.975),quantile(full_inc_ME,0.025),'b','None',1,0.2); hold on;
pl5=plot(1:tf, mean(full_inc_LO),'linewidth',2.5,'color',[0.90,0.35,0.11]);hold on;
% errorbar(find(mean(full_inc_LO)==max(mean(full_inc_LO))), max(mean(full_inc_LO)), ...
%     max(mean(full_inc_LO))-max(quantile(full_inc_LO,0.025)), max(quantile(full_inc_LO,0.975))-max(mean(full_inc_LO)), 'color',[0.90,0.35,0.11])
%jbfill(1:tf,quantile(full_inc_LO,0.975),quantile(full_inc_LO,0.025),[0.90,0.35,0.11],'None',1,0.2); hold on;
pl6=plot(1:tf, mean(full_inc_HI),'linewidth',2.5,'color',[0.49,0.18,0.56]);hold on;
% errorbar(find(mean(full_inc_HI)==max(mean(full_inc_HI))), max(mean(full_inc_HI)), ...
%     max(mean(full_inc_HI))-max(quantile(full_inc_HI,0.025)), max(quantile(full_inc_HI,0.975))-max(mean(full_inc_HI)), 'color',[0.49,0.18,0.56])
%jbfill(1:tf,quantile(full_inc_HI,0.975),quantile(full_inc_HI,0.025),[0.49,0.18,0.56],'None',1,0.2); hold on;

pl2=plot(1:tf,mean(tinc4),'linewidth',2.5,'color','r');
% errorbar(find(mean(tinc4)==max(mean(tinc4))), max(mean(tinc4)), ...
%     max(mean(tinc4))-max(quantile(tinc4,0.025)), max(quantile(tinc4,0.975))-max(mean(tinc4)), 'color','r')
%jbfill(1:tf,quantile(tinc4,0.975),quantile(tinc4,0.025),'r','None',1,0.2); hold on;
pl3=plot(1:tf,mean(tinc5),'linewidth',2.5,'color',[0.47,0.67,0.19]);
% errorbar(find(mean(tinc5)==max(mean(tinc5))), max(mean(tinc5)), ...
%     max(mean(tinc5))-max(quantile(tinc5,0.025)), max(quantile(tinc5,0.975))-max(mean(tinc5)), 'color',[0.47,0.67,0.19])
%jbfill(1:tf,quantile(tinc5,0.975),quantile(tinc5,0.025),[0.47,0.67,0.19],'None',1,0.2);hold on;
pl4=plot(1:tf,mean(tinc0),'linewidth',2.5,'color','black');
% errorbar(find(mean(tinc0)==max(mean(tinc0))), max(mean(tinc0)), ...
%     max(mean(tinc0))-max(quantile(tinc0,0.025)), max(quantile(tinc0,0.975))-max(mean(tinc0)), 'color','black')
jbfill(1:tf,quantile(tinc0,0.975),quantile(tinc0,0.025),'black','None',1,0.2);
hold off; box on;
legend([pl4 pl2 pl3 pl1 pl5 pl6],{'Baseline (no intervention)','LFA testing','LFA+PCR testing',...
    'Dynamic testing (medium threshold)', 'Dynamic testing (low threshold)','Dynamic testing (high threshold)'});
% legend([pl4 pl2 pl3 pl1 pl5 pl6],{'Baseline (no intervention)','LFA testing','LFA+PCR testing',...
%     'Dynamic testing (medium threshold, RDT)', 'Dynamic testing (low threshold, RDT)','Dynamic testing (high threshold, RDT)'});
set(gca,'fontsize',18);
xlabel('Days'); ylabel('Daily symptomatic incidence (thousands)')
%ylabel({'Daily symptomatic'; 'incidence (thousands)'});


figure;
% -- plot PCR and LFA consumption over the intervention period -- %
subplot(1,2,1); 
% xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
% xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(0:tf,mean(pcr_cons_ME),'linewidth',2.5); hold on;
pl4=plot(0:tf,mean(pcr_cons_LO),'linewidth',2.5,'color',[0.90,0.35,0.11]);
pl5=plot(0:tf,mean(pcr_cons_HI),'linewidth',2.5,'color',[0.49,0.18,0.56]);
pl2=plot(0:tf,mean(pcr5),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
% pl3=plot(0:tf,mean(pcr0),'linewidth',2.5,'color','black'); 
hold off; title('PCR consumption');
legend([pl2 pl1 pl4 pl5],{'LFA+PCR strategy',...
    'Dynamic testing (medium threshold)', 'Dynamic testing (low threshold)','Dynamic testing (high threshold)'}); hold off;
set(gca,'fontsize',16); ylabel('Cumulative number of tests performed (millions)');
xlabel('Days'); box on;

% subplot(1,2,2);
% % xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
% % xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
% pl1=plot(0:tf,mean(lfa_cons_ME),'linewidth',2.5); hold on;
% pl4=plot(0:tf,mean(lfa_cons_LO),'linewidth',2.5,'color',[0.90,0.35,0.11]);
% pl5=plot(0:tf,mean(lfa_cons_HI),'linewidth',2.5,'color',[0.49,0.18,0.56]);
% pl2=plot(0:tf,mean(lfa4),'linewidth',2.5);
% pl3=plot(0:tf,mean(lfa5),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
% title('LFA consumption'); 
% hold off; ylim([0 max(max(lfa4))]);
% legend([pl3 pl2 pl1 pl4 pl5],{'Baseline (no intervention)','LFA+PCR strategy',...
%     'Dynamic testing (medium threshold)', 'Dynamic testing (low threshold)','Dynamic testing (high threshold)'}); hold off;
% set(gca,'fontsize',16); ylabel({'Cumulative number of','tests performed (millions)'});
% xlabel('Days'); box on;

% -- False positive quars
subplot(1,2,2);
% xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
% xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(0:tf,mean(fp_quar_ME)./sum(prm.N),'linewidth',2.5); hold on;
pl4=plot(0:tf,mean(fp_quar_LO)./sum(prm.N),'linewidth',2.5,'color',[0.90,0.35,0.11]);
pl5=plot(0:tf,mean(fp_quar_HI)./sum(prm.N),'linewidth',2.5,'color',[0.49,0.18,0.56]);
pl2=plot(0:tf,mean(fpqur4)./sum(prm.N),'linewidth',2.5,'color','r');
pl3=plot(0:tf,mean(fpqur5)./sum(prm.N),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
title('Unnecessary isolations'); 
hold off; %ylim([0 max(max(lfa4))]);
legend([pl2 pl3 pl1 pl4 pl5],{'LFA strategy','LFA+PCR strategy',...
    'Dynamic testing (medium threshold)', 'Dynamic testing (low threshold)','Dynamic testing (high threshold)'}); hold off;
set(gca,'fontsize',16); ylabel('Unnecessary isolations (relative to population size)');
xlabel('Days'); box on;

% -- Boxplots
figure; 
subplot(1,2,1);title('PCR consumption');
boxplot([ones(150,1)*nan, ca_ME(:,3)./pcr5(:,end),ca_ME(:,1)./pcr_cons_ME(:,end),ca_LO(:,1)./pcr_cons_LO(:,end),ca_HI(:,1)./pcr_cons_HI(:,end)],...
    'Labels',{'LFA strategy','LFA+PCR strategy','DY (MED)','DY (LOW)','DY (HI)'},...
    'Notch','on');
text(.9,0.05, "NA","FontSize",18)
title('Cases averted per PCR test used');
box on; set(gca,'fontsize',16); 
subplot(1,2,2); title('Unnecessary isolations');
boxplot([ca_ME(:,2)./fpqur4(:,end),ca_ME(:,3)./fpqur5(:,end),ca_ME(:,1)./fp_quar_ME(:,end),ca_LO(:,1)./fp_quar_LO(:,end),ca_HI(:,1)./fp_quar_HI(:,end)],...
    'Labels',{'LFA strategy','LFA+PCR strategy','DY (MED)','DY (LOW)','DY (HI)'},...
    'Notch','on');
title('Cases averted per unnecessary isolation');
box on;
set(gca,'fontsize',16); 

% figure; 
% title('Unnecessary isolations');
% boxplot([ca_ME(:,2)./fpqur4(:,end),ca_ME(:,1)./fp_quar_ME(:,end),ca_LO(:,1)./fp_quar_LO(:,end),ca_HI(:,1)./fp_quar_HI(:,end)],...
%     'Labels',{'LFA strategy','DY (MED)','DY (LOW)','DY (HI)'},...
%     'Notch','on');
% title('Cases averted per unnecessary isolation');
% box on;
% set(gca,'fontsize',16); 
