%% Calibrate non-covid population to coincide with second wave, to give 30% positivity rate in the clinics at peak

clear all; load Model_setup_June.mat;

R0 = 2.5;
r.beta = r.beta*R0;

fac.lo = 1/2;
fac.md = 1;
fac.hi = 2;
opt = 'md'; % 'lo', 'md','hi'

% --- Control parameters we may want to stratify by -----------------------
prm.pNCS    = 0.09;
r.careseek  = 0.1;
r.Zcareseek = 0.02;
p_seropos   = 0.3; %% second wave
nrdt        = 0;%1.2466e+07*15*fac.(opt);
% nrdt       = 1e11;

% prm.PCR_capacity = inf;

% Set 'testing' to 1 if you don't want to do full uncertainty (takes ~2hrs 
% for 100 samples), and only want to look at some example results with a 
% single parameter set (takes 1-2 mins). Set to 0 otherwise
testing = 1;



tf = 600;
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


isam=1;
[p,r,prm] = alloc_parameters3(xsam(isam,:), psto, rsto, xi, i, s, gps, prmsto, R0);

init = zeros(1,i.nx); seed = 10;
init(intersect(intersect(s.S,s.q0),s.null)) = prm.N*(1-p_seropos);
init(intersect(intersect(s.R,s.q0),s.null)) = prm.N*p_seropos;
init(i.I.null.q0.ad) = seed; init(i.S.null.q0.ad) = init(i.S.null.q0.ad) - seed;
init(i.Z.null.q0)        = sum(prm.N)*prm.pNCS;


% --- Baseline
p0 = p; r0 = r;
p0.LFA_cs = 0;
r0.screen = zeros(3,2);
M0 = make_model_iterative(p0, r0, i, s, gps, prm);
geq = @(t,in) goveqs_basis(t, in, M0, i, s, r0, prm, agg, sel);
[~,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));
cinc0 = sum(soln0(end,i.aux.inc),2);
tinc  = sum(diff(soln0(:,i.aux.inc),1),2);

soln0(end,i.aux.pcr)
soln0(end,i.aux.pcr(1))/sum(soln0(end,i.aux.pcr))

% figure; plot(tinc,'linewidth',2.5);
% set(gca,'fontsize',16);
% 
% return;
% Check COVID prevalence in clinics
figure; subplot(1,2,1);
plot(sum(soln0(:,intersect(s.I,s.pcra)), 2),'linewidth',2.5);ylabel('COVID-19 cases self-presenting for care');
xlabel('Days'); set(gca,'fontsize',16); xlim([0 tf]);
subplot(1,2,2); 
plot(sum(soln0(:,intersect(s.I,s.pcra)), 2)./(sum(soln0(:,intersect(s.I,s.pcra)), 2)+ sum(soln0(:,intersect(s.Z,s.pcra)), 2)),'linewidth',2.5)
ylabel('Prevalence of COVID-19 in clinics');xlabel('Days');xlim([0 tf]);
set(gca,'fontsize',16);



