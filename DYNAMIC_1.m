% givencov: Simulate epi outcomes under given coverage rates in clinic and
% community
set(0,'DefaultFigureWindowStyle','docked')
clear all; load Model_setup_Feb22.mat;

R0 = 2.5;
r.beta = r.beta*R0;
% cost = [50 5];

prm.pNCS            = 0.09;
r.careseek          = 0.1;
r.Zcareseek         = 0.1;
p_seropos           = 0.3;
prm.PCR_capacity    = 1e10; % Inf
tf = 600;

r_screen   = 6e7/sum(prm.N)/tf;
p.LFA_cs   = 1;

testing = 0;

% r.Dx(1)= 1/30;    % 1 week - 1 month delay before being receiving PCR confirmation
p.spec(2) = 0.99; % increased specificity of LFAs, now equal to PCR spec
r.hold    = 1e4;

psto = p; rsto = r; prmsto = prm;

% --- Get random samples within uncertainty bounds ------------------------
if testing
    % Just take a single, example parameter set
    nsam = 1;
    xsam = mean(prm.bounds,1);
else
    % Do full uncertainty sampling
    nsam = 150;
    xsam = repmat(prm.bounds(1,:),nsam,1) + lhsdesign(nsam, xi.nx).*repmat(diff(prm.bounds,1),nsam,1);
end

% threshold = 0.50; %%% ON/OFF

tic
mk = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end
    
    [p,r,prm] = alloc_parameters3(xsam(ii,:), psto, rsto, xi, i, s, gps, prmsto, R0);
    
    % Setting up initial conditions
    p_seropos = 0.3;
    init = zeros(1,i.nx); seed = 10;
    init(intersect(intersect(s.S,s.q0),s.null)) = prm.N*(1-p_seropos);
    init(intersect(intersect(s.R,s.q0),s.null)) = prm.N*p_seropos;
    %init(i.I.null.q0.ad.mld) = seed; init(i.S.null.q0.ad) = init(i.S.null.q0.ad) - seed;
    init(i.I.null.q0.ad) = seed; init(i.S.null.q0.ad) = init(i.S.null.q0.ad) - seed;
    init(i.Z.null.q0)    = sum(prm.N)*prm.pNCS;
        
    % --- Infection is low, we use LFAs+PCR for confirmation ----------------------
    p0 = p; r0 = r;
    p0.LFA_cs  = 1;
    r0.Dx(1)   = 1/3;    % 3 day delay; 1 week - 1 month delay before being receiving PCR confirmation
    r0.screen  = [zeros(3,1), r_screen*ones(3,1)];
    M0 = make_model_iterative_dynamic(p0, r0, i, s, gps, prm);
    geq = @(t,in) goveqs_basis(t, in, M0, i, s, r0, prm, agg, sel);
    [~,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
%     swi = soln0(:,i.aux.pcr(1))./sum(soln0(:,i.aux.pcr),2); % chosen switch: proportion of true +ve tests

    swi = sum(diff(soln0(:,i.aux.qur),1),2)./sum(diff(soln0(:,i.aux.pcr),1),2); % chosen switch: test positivity - PCR (positive / all tests)
%     swi = sum(diff(soln0(:,i.aux.pcr),1),2)./sum(diff(soln0(:,i.aux.rdt),1),2); % chosen switch: test positivity - RDT (positive / all tests)
    
    % -- save test positivity under baseline scenario for plotting purposes later
    pcr_pos_rate(ii,:) = sum(diff(soln0(:,i.aux.qur),1),2)./sum(diff(soln0(:,i.aux.pcr),1),2);
    rdt_pos_rate(ii,:) = sum(diff(soln0(:,i.aux.pcr),1),2)./sum(diff(soln0(:,i.aux.rdt),1),2);

    tinc = sum(diff(soln0(:,i.aux.inc),1),2);
    fii  = find(swi>max(swi)*.1);% chosen threshold %%% max(swi)*.5 / max(swi)*.9 / max(swi)*.1
    t1   = fii(1);   
    fpqur1  = sum(soln0(:,i.aux.fpqur),2); % @t1 + @t2 +
    pcr1    = sum(soln0(:,i.aux.pcr),2); % PCR CONSUMPTION
    lfa1    = sum(soln0(:,i.aux.rdt),2); % LFA CONSUMPTION
    
    % --- Epidemic has taken off, we drop PCR confirmation --------------------------------------------    
    p1 = p0; r1 = r0;
    r1.Dx(1)   = r.Dx(1);
    p1.sens(2) = p.sens(2); 
    p1.spec(2) = p.spec(2);
    M1 = make_model_iterative(p1, r1, i, s, gps, prm); %% Use a diff make_model with different transitions through testing compartments
    geq = @(t,in) goveqs_basis(t, in, M1, i, s, r1, prm, agg, sel);
    [~,soln1] = ode15s(geq, [t1:1:tf], soln0(t1,:), odeset('Nonnegative',1:i.nx));
    
    tinc2 = sum(diff(soln1(:,i.aux.inc),1),2); %%%
%     swi_off = sum(soln1(:,i.aux.fpqur)./sum(prm.N),2); % switch off variable

    % - switch off when test positivity is above threshold again, but here we look at RDT positivity because there's no PCR   
    swi = sum(diff(soln1(:,i.aux.qur),1),2)./sum(diff(soln1(:,i.aux.rdt),1),2);
%     fii = find(swi_off>threshold);% chosen threshold %%%
    fii = find(swi<max(swi)*.3);% chosen threshold %%% max(swi)*.5 / max(swi(3:end))*.9 / max(swi)*.3
    fii = fii(fii>=find(tinc2>=max(tinc2))); %% How to fix this?? Atm i'm manually restricting the switch back on to happen after epi peak %%
    t2  = t1+fii(1);    
    fpqur2  = sum(soln1(:,i.aux.fpqur),2); % @t1 + @t2 + 
    pcr2    = sum(soln1(:,i.aux.pcr),2); % PCR CONSUMPTION
    lfa2    = sum(soln1(:,i.aux.rdt),2); % LFA CONSUMPTION

    % --- Epi down again, back to PCR confirmation -------------------------------------------
    p2 = p0; r2 = r0;
    M2 = make_model_iterative_dynamic(p2, r2, i, s, gps, prm);
    geq = @(t,in) goveqs_basis(t, in, M2, i, s, r2, prm, agg, sel);
    [~,soln2] = ode15s(geq, [t2:1:tf], soln1(fii(1),:), odeset('Nonnegative',1:i.nx));
    tinc3   = sum(diff(soln2(:,i.aux.inc),1),2);
    fpqur3  = sum(soln2(:,i.aux.fpqur),2); % @t1 + @t2 + 
    pcr3    = sum(soln2(:,i.aux.pcr),2); % PCR CONSUMPTION
    lfa3    = sum(soln2(:,i.aux.rdt),2); % LFA CONSUMPTION
    cinc_dy = sum(soln2(end,i.aux.inc));

    % --- Bring all solutions together ------------------------------------
%     full_sol = [soln0(1:t1,:), soln1(t1+1:t2,:), soln2(t2+1:tf)];
    full_inc(ii,:) = [tinc(1:t1)', tinc2(1:fii(1))', tinc3'];
    fp_quar(ii,:)  = [fpqur1(1:t1)', fpqur2(1:fii)', fpqur3']; 
    pcr_cons(ii,:) = [pcr1(1:t1)', pcr2(1:fii)', pcr3'];
    lfa_cons(ii,:) = [lfa1(1:t1)', lfa2(1:fii)', lfa3'];
    
    
    % --- Other scenarios
    % --- LFAs only in community and clinics -------------------------------------------
    p3 = p; r3 = r;
    r3.screen = [zeros(3,1), r_screen*ones(3,1)];
    p3.LFA_cs = 1;
    M3 = make_model_iterative(p3, r3, i, s, gps, prm);
    geq = @(t,in) goveqs_basis(t, in, M3, i, s, r3, prm, agg, sel);
    [~,soln3] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
    
    tinc4(ii,:) = sum(diff(soln3(:,i.aux.inc),1),2);
    fpqur4(ii,:)= sum(soln3(:,i.aux.fpqur),2); % 
    lfa4(ii,:)  = sum(soln3(:,i.aux.rdt),2); % LFA CONSUMPTION
    cinc4       = sum(soln3(end,i.aux.inc));

    % --- LFAs+PCR in community and clinics throughout the timeline ---------------------
    p4 = p; r4 = r;
    r4.Dx(1) = r2.Dx(1); % delay of PCR confirmation
    r4.screen = [zeros(3,1), r_screen*ones(3,1)];
    p4.LFA_cs = 1;
    M4 = make_model_iterative_dynamic(p4, r4, i, s, gps, prm);
    geq = @(t,in) goveqs_basis(t, in, M4, i, s, r4, prm, agg, sel);
    [~,soln4] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
    
    tinc5(ii,:) = sum(diff(soln4(:,i.aux.inc),1),2);
    fpqur5(ii,:)= sum(soln4(:,i.aux.fpqur),2); % 
    pcr5(ii,:)  = sum(soln4(:,i.aux.pcr),2); % PCR CONSUMPTION
    lfa5(ii,:)  = sum(soln4(:,i.aux.rdt),2); % LFA CONSUMPTION
    cinc5       = sum(soln4(end,i.aux.inc));
    
    % --- Baseline, no intervention with LFAs
    p0 = p; r0 = r;
    p0.LFA_cs = 0;
    r_screen  = 0;%6e7/sum(prm.N)/tf; 
    r0.screen = [zeros(3,1), r_screen*ones(3,1)];
    prm1 = prm; prm1.PCR_capacity = 1e5;
    M_base = make_model_iterative(p0, r0, i, s, gps, prm1);
    geq = @(t,in) goveqs_basis(t, in, M_base, i, s, r0, prm1, agg, sel);
    [~,soln_base] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
    tinc0(ii,:)= sum(diff(soln_base(:,i.aux.inc),1),2);
    fpqur0     = sum(soln_base(t1,i.aux.fpqur)); % @t1 + @t2 + 
    cinc0      = sum(soln_base(end,i.aux.inc));
    pcr0(ii,:) = sum(soln_base(:,i.aux.pcr),2); % PCR CONSUMPTION

    
    pca(ii,1) = 1 - cinc_dy/cinc0;
    pca(ii,2) = 1 - cinc4/cinc0;
    pca(ii,3) = 1 - cinc5/cinc0;    
    
    ca(ii,1)  = cinc0 - cinc_dy;
    ca(ii,2)  = cinc0 - cinc4;
    ca(ii,3)  = cinc0 - cinc5;
    
end
toc

fprintf('\n');

% figure; subplot(1,2,1);
% plot(1:tf,mean(pcr_pos_rate)*100,'linewidth',2.5); hold on;
% yline([max(mean(pcr_pos_rate)) max(mean(pcr_pos_rate))]*90,'--','linewidth',2,'color','red',alpha=0.2);
% yline([max(mean(pcr_pos_rate)) max(mean(pcr_pos_rate))]*50,'--','linewidth',2,'color','magenta',alpha=0.2);
% yline([max(mean(pcr_pos_rate)) max(mean(pcr_pos_rate))]*10,'--','linewidth',2,'color','cyan',alpha=0.2);
% hold off;
% xlabel('Days'); box on;
% ylabel('Test positivity rate (PCR)');
% set(gca,'fontsize',16);
% 
% subplot(1,2,2); plot(1:tf,mean(rdt_pos_rate)*100,'linewidth',2.5);
% hold on;
% yline([max(mean(rdt_pos_rate)) max(mean(rdt_pos_rate))]*90,'--','linewidth',2,'color','red',alpha=0.2);
% yline([max(mean(rdt_pos_rate)) max(mean(rdt_pos_rate))]*50,'--','linewidth',2,'color','magenta',alpha=0.2);
% yline([max(mean(rdt_pos_rate)) max(mean(rdt_pos_rate))]*30,'--','linewidth',2,'color','cyan',alpha=0.2);
% hold off;
% xlabel('Days'); box on;
% ylabel('Test positivity rate (RDT)');
% set(gca,'fontsize',16);

save('dynamic_results_LO_221121.mat'); % ME, LO, HI for the different thresholds % RDT for RDT test pos analysis
return;

figure; subplot(2,3,[1 2 3]);
xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(1:tf, mean(full_inc),'linewidth',2.5);hold on;
jbfill(1:tf,quantile(full_inc,0.975),quantile(full_inc,0.025),'b','None',1,0.2); hold on;
% line([1,600],[threshold,threshold],'linestyle','--')
pl2=plot(1:tf,mean(tinc4),'linewidth',2.5);
jbfill(1:tf,quantile(tinc4,0.975),quantile(tinc4,0.025),'r','None',1,0.2); hold on;
pl3=plot(1:tf,mean(tinc5),'linewidth',2.5,'color',[0.47,0.67,0.19]);
jbfill(1:tf,quantile(tinc5,0.975),quantile(tinc5,0.025),[0.47,0.67,0.19],'None',1,0.2);hold on;
pl4=plot(1:tf,mean(tinc0),'linewidth',2.5,'color','black');
jbfill(1:tf,quantile(tinc0,0.975),quantile(tinc0,0.025),'black','None',1,0.2);
hold off; box on;
legend([pl1 pl2 pl3 pl4],{'Dynamic testing','LFA testing','LFA+PCR testing','Baseline (no intervention)'});
set(gca,'fontsize',16);
xlabel('Days'); ylabel({'Daily symptomatic'; 'incidence (thousands)'});
% return;

% -- plot PCR and LFA consumption over the intervention period -- %
subplot(2,3,4); 
xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(0:tf,mean(pcr_cons),'linewidth',2.5);
pl2=plot(0:tf,mean(pcr5),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
pl3=plot(0:tf,mean(pcr0),'linewidth',2.5,'color','black'); 
hold off; title('PCR consumption');
legend([pl1 pl2 pl3],{'Dynamic strategy','LFA+PCR strategy','Baseline (no intervention)'}); hold off;
set(gca,'fontsize',16); ylabel({'Cumulative number of','tests performed (thousands)'});
xlabel('Days'); box on;

subplot(2,3,5);
xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(0:tf,mean(lfa_cons),'linewidth',2.5);
pl2=plot(0:tf,mean(lfa4),'linewidth',2.5);
pl3=plot(0:tf,mean(lfa5),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
title('LFA consumption'); 
hold off; ylim([0 max(max(lfa4))]);
legend([pl1 pl2 pl3],{'Dynamic strategy','LFA strategy','LFA+PCR strategy'});
set(gca,'fontsize',16); ylabel({'Cumulative number of','tests performed (millions)'});
xlabel('Days'); box on;

% -- False positive quars
subplot(2,3,6);
xline([t1 t1],'--','linewidth',2,'color','red',alpha=0.2);hold on;
xline([t2 t2],'--','linewidth',2,'color','red',alpha=0.2);
pl1=plot(0:tf,mean(fp_quar)./sum(prm.N),'linewidth',2.5);
pl2=plot(0:tf,mean(fpqur4)./sum(prm.N),'linewidth',2.5);
pl3=plot(0:tf,mean(fpqur5)./sum(prm.N),'linewidth',2.5,'color',[0.47,0.67,0.19]); 
title('Unnecessary isolations'); 
hold off; %ylim([0 max(max(lfa4))]);
legend([pl1 pl2 pl3],{'Dynamic strategy','LFA strategy','LFA+PCR strategy'});
set(gca,'fontsize',16); ylabel({'Unnecessary isolations'; '(relative to population size)'});
xlabel('Days'); box on;


return;
%% Varying the incidence threshold of the dynamic model
% clear all;
thresh = [.25:.05:1]; %

tic
mk = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end

    [p,r,prm] = alloc_parameters3(xsam(ii,:), psto, rsto, xi, i, s, gps, prmsto, R0);
    
    % Setting up initial conditions
    p_seropos = 0.3;
    init = zeros(1,i.nx); seed = 10;
    init(intersect(intersect(s.S,s.q0),s.null)) = prm.N*(1-p_seropos);
    init(intersect(intersect(s.R,s.q0),s.null)) = prm.N*p_seropos;
    %init(i.I.null.q0.ad.mld) = seed; init(i.S.null.q0.ad) = init(i.S.null.q0.ad) - seed;
    init(i.I.null.q0.ad) = seed; init(i.S.null.q0.ad) = init(i.S.null.q0.ad) - seed;
    init(i.Z.null.q0)        = sum(prm.N)*prm.pNCS;
    
    for th = 1:length(thresh)
            threshold = thresh(th);

            % --- Infection is low, we use LFAs+PCR for confirmation ----------------------
            p0 = p; r0 = r;
            p0.LFA_cs  = 1;
            r0.Dx(1)   = 1/7;    % 1 week - 1 month delay before being receiving PCR confirmation
            r0.screen  = [zeros(3,1), r_screen*ones(3,1)];
            M0 = make_model_iterative_dynamic(p0, r0, i, s, gps, prm);
            geq = @(t,in) goveqs_basis(t, in, M0, i, s, r0, prm, agg, sel);
            [~,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
        %     swi = soln0(:,i.aux.pcr(1))./sum(soln0(:,i.aux.pcr),2); % chosen switch: proportion of true +ve tests
        
            swi = sum(diff(soln0(:,i.aux.qur),1),2)./sum(diff(soln0(:,i.aux.pcr),1),2); % chosen switch: test positivity - PCR (positive / all tests)
        %     swi = sum(diff(soln0(:,i.aux.pcr),1),2)./sum(diff(soln0(:,i.aux.rdt),1),2); % chosen switch: test positivity - RDT (positive / all tests)
            
            % -- save test positivity under baseline scenario for plotting purposes later
            pcr_pos_rate(ii,:) = sum(diff(soln0(:,i.aux.qur),1),2)./sum(diff(soln0(:,i.aux.pcr),1),2);
            rdt_pos_rate(ii,:) = sum(diff(soln0(:,i.aux.pcr),1),2)./sum(diff(soln0(:,i.aux.rdt),1),2);
        
            tinc = sum(diff(soln0(:,i.aux.inc),1),2);
            fii  = find(swi>max(swi)*threshold);% chosen threshold %%% max(swi)*.5 / max(swi)*.9 / max(swi)*.1
            
            if isempty(fii) % only when threshold == 1
                %t1      = 600; % never switch off PCR confirmation
                cinc_dy        = sum(soln0(end,i.aux.inc)); 
                fp_quar(ii,:)  = sum(soln0(:,i.aux.fpqur),2);
                % Save outputs
                pca_dyn(ii,th)  = 1 - cinc_dy/cinc0;
                inc_dy(th,:)    = full_inc(ii,:);
                
                % cases averted per unnecessary isolation, over different switch thresholds
                ca_av(ii,th)    = (cinc0 - cinc_dy)/fp_quar(ii,end);
            else
                t1      = fii(1);   
                fpqur1  = sum(soln0(:,i.aux.fpqur),2); % @t1 + @t2 +
                
                % --- Epidemic has taken off, we drop PCR confirmation --------------------------------------------    
                p1 = p0; r1 = r0;
                r1.Dx(1)   = r.Dx(1);
                p1.sens(2) = p.sens(2); 
                p1.spec(2) = p.spec(2);
                M1 = make_model_iterative(p1, r1, i, s, gps, prm); %% Use a diff make_model with different transitions through testing compartments
                geq = @(t,in) goveqs_basis(t, in, M1, i, s, r1, prm, agg, sel);
                [~,soln1] = ode15s(geq, [t1:1:tf], soln0(t1,:), odeset('Nonnegative',1:i.nx));
                
                tinc2 = sum(diff(soln1(:,i.aux.inc),1),2); %%%
                % - switch off when test positivity is above threshold again, but here we look at RDT positivity because there's no PCR   
                swi = sum(diff(soln1(:,i.aux.qur),1),2)./sum(diff(soln1(:,i.aux.rdt),1),2);
                fii = find(swi<max(swi)*threshold);% chosen threshold 
                fii = fii(fii>=find(tinc2>=max(tinc2))); 
                t2  = t1+fii(1);    
                fpqur2  = sum(soln1(:,i.aux.fpqur),2); % @t1 + @t2 + 
            
                % --- Epi down again, back to PCR confirmation -------------------------------------------
                p2 = p0; r2 = r0;
                M2 = make_model_iterative_dynamic(p2, r2, i, s, gps, prm);
                geq = @(t,in) goveqs_basis(t, in, M2, i, s, r2, prm, agg, sel);
                [~,soln2] = ode15s(geq, [t2:1:tf], soln1(fii(1),:), odeset('Nonnegative',1:i.nx));
                tinc3   = sum(diff(soln2(:,i.aux.inc),1),2);
                fpqur3  = sum(soln2(:,i.aux.fpqur),2); % @t1 + @t2 + 
                cinc_dy = sum(soln2(end,i.aux.inc));
            
                % --- Bring all solutions together ------------------------------------
            %     full_sol = [soln0(1:t1,:), soln1(t1+1:t2,:), soln2(t2+1:tf)];
                full_inc(ii,:) = [tinc(1:t1)', tinc2(1:fii(1))', tinc3'];
                fp_quar(ii,:)  = [fpqur1(1:t1)', fpqur2(1:fii)', fpqur3'];                 
            
            
                % --- Baseline, no intervention with LFAs
                p_b = p; r_b = r;
                p_b.LFA_cs = 0;
                r_screen   = 0;%6e7/sum(prm.N)/tf; 
                r_b.screen = [zeros(3,1), r_screen*ones(3,1)];
                prm1 = prm; prm1.PCR_capacity = 1e5;
                M_base = make_model_iterative(p_b, r_b, i, s, gps, prm1);
                geq = @(t,in) goveqs_basis(t, in, M_base, i, s, r_b, prm1, agg, sel);
                [~,soln_base] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
                tinc0(ii,:)= sum(diff(soln_base(:,i.aux.inc),1),2);
                fpqur0     = sum(soln_base(t1,i.aux.fpqur)); % @t1 + @t2 + 
                cinc0      = sum(soln_base(end,i.aux.inc));
    
    
                % Save outputs
                pca_dyn(ii,th)  = 1 - cinc_dy/cinc0;
                inc_dy(th,:)    = full_inc(ii,:);
                
                % cases averted per unnecessary isolation, over different switch thresholds
                ca_av(ii,th)    = (cinc0 - cinc_dy)/fp_quar(ii,end);
            end
    end
end

%save('switch_threshold.mat');
figure; 
plot(thresh,mean(ca_av),'linewidth',2.5,'Marker','o'); hold on;
jbfill(thresh,quantile(ca_av,0.975),quantile(ca_av,0.025),'b','None',1,0.2); hold off;
xlabel('Switch threshold (% of the maximum test positivity rate)'); 
ylabel('Cases averted per unnecessary isolation');
set(gca,'fontsize',14);
xlim([0.25,1]); ylim([0, max(max(ca_av))*1.1]);

return;

% figure; subplot(1,2,1);
% plot(thresh,mean(pca_dyn)*100,'linewidth',2.5); %title()
% hold on;
% jbfill(thresh,quantile(pca_dyn,0.975)*100,quantile(pca_dyn,0.025)*100,'b','None',1,0.2); hold off;
% ylabel('Percent cases averted'); ylim([0 max(max(pca_dyn))*100]);
% xlabel('PCR testing delay (days)'); xlim([thresh(1) thresh(end)])
% set(gca,'fontsize',14);
% 
% subplot(1,2,2);
% plot(1:tf,inc_dy(1,:),'linewidth',2.5); hold on;
% plot(1:tf,inc_dy(2,:),'linewidth',2.5);
% plot(1:tf,inc_dy(3,:),'linewidth',2.5);
% plot(1:tf,inc_dy(4,:),'linewidth',2.5);
% xlabel('Days'); 
% ylabel('Daily symptomatic incidence (thousands)');
% legend({'1 week','2 weeks','3 weeks','4 weeks'});
% set(gca,'fontsize',14);
