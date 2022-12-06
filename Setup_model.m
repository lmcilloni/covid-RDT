% v2: Limiting to just one city (Delhi), to help simulation time

clear all; 

gps = [];

gps.age      = {'ch','ad','el'};
gps.quar     = {'q0','q1'};
gps.test     = {'null','pcra','pcrb','rdt'};

states1  = {'S','E','A','P','R'};                                          % Age, quar status, testing status (PCR/RDT)
states2  = {'I'};                                                          % Age, quar status, testing status
states3  = {'D'};                                                          % Age only
states4  = {'Z'};                                                          % Quar status, testing status

% [i, s, d, lim] = get_addresses({states1, gps.age, gps.quar, gps.test}, [], [], [], 0);
% [i, s, d, lim] = get_addresses({states2, gps.age, gps.quar, gps.test, gps.severity}, i, s, d, lim);
% [i, s, d, lim] = get_addresses({states3, gps.age}, i, s, d, lim);
% s.main = 1:lim;
% [i, s, d, lim] = get_addresses({states4, gps.quar, gps.test}, i, s, d, lim);
% d = char(d);

[i, s, d, lim] = get_addresses({states1, gps.test, gps.quar, gps.age}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.test, gps.quar, gps.age}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states3, gps.age}, i, s, d, lim);
s.main = 1:lim;
[i, s, d, lim] = get_addresses({states4, gps.test, gps.quar}, i, s, d, lim);
d = char(d);

% This step is to create numerical indices for each of the state variables,
% which makes it easy to set up the governing equations when the number of
% state variables becomes large. This approach depends heavily on Matlab 
% 'structures', which are the same as 'lists' in R. 
% In the above lines:
% - 's' contains the numerical indices of different groups of states. For
% example, s.E gives you all stratifications of E, and s.mld give you all
% states relating to mild illness (see pdf for list of states)
% - 'i' gives the numerical indices of specific states. For example, type
% i.S.el to find that the index of susceptible, elderly people is 3.
% - 'd' is a full list of the state variables, in numerical order

s.infectious = intersect([s.A,s.P,s.I],s.q0);                              % Indices of all infectious compartments
s.prevalent  = [s.I];                                                 % Indices of all symptomatic prevalence compartments 
s.Sq0        = intersect(s.S, s.q0);                                       % Indices of all susceptible who are not in quarantine
s.Eq0        = intersect(s.E, s.q0);                                       % Indices of all exposed who are not in quarantine

% % s.prevalent  = [s.I, s.Dx1, s.Dx2];                                        
% % s.S1         = intersect(s.S,s.q0);                                        % Susceptibles liable to infection
% % s.E1         = intersect(s.E,s.q0);                                        % When such susceptibles are infected, the state they enter

% --- Include the auxiliaries ---------------------------------------------

numcoms = length(gps.age);
names   = {  'inc',  'pcr', 'rdt_comm', 'rdt_clin', 'fpqur','qur'};
lgths   = [numcoms,      2,          2,          2,       2,    2];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.aux.rdt = [i.aux.rdt_clin, i.aux.rdt_comm];
i.nx = lim;

% The above are additional variables that we include in the governing
% equations, to help count 'incidence-like' terms, e.g. the daily
% symptomatic incidence, and the daily admissions to hospital. They're not
% state variables, instead we refer to them as 'auxiliaries', with the
% label 'aux' incorporated in their numerical indices. When we include 
% counts of the numbers of tests being performed per day, they'll be 
% included here as more 'incidence-like' terms. 


% --- Counting symptomatic incidence, by age
mat = zeros(i.nstates);
mat(s.I,s.P)   = 1;
mat(s.q1,s.q0) = 0; mat(s.q0,s.q1) = 0;
sel.inc = sparse(mat - diag(diag(mat)));

mat = zeros(numcoms,i.nstates);
for ia = 1:length(gps.age)
    mat(ia, intersect(s.I,s.(gps.age{ia}))) = 1;
end
agg.inc = sparse(mat);


% % --- Counting hospitalisations
% mat = zeros(i.nstates);
% mat(s.H,:) = 1;
% sel.hosp = sparse(mat - diag(diag(mat)));
% 
% mat = zeros(numcoms,i.nstates); row = 1;
% for ia = 1:length(gps.age)
%     mat(ia, intersect(s.H,s.(gps.age{ia}))) = 1;
% end
% agg.hosp = sparse(mat);


% --- Counting PCR tests being used
mat = zeros(i.nstates);
mat(s.pcrb,s.pcra)  = 1;
sel.pcr = sparse(mat - diag(diag(mat)));

mat = zeros(2,i.nstates);
mat(1,intersect(s.main,s.pcrb)) = 1;
mat(2,intersect(s.Z,s.pcrb))    = 1;
agg.pcr = sparse(mat);


% --- Counting RDT tests being used
mat = zeros(i.nstates);
mat(s.rdt,s.null)  = 1;
sel.rdt = sparse(mat - diag(diag(mat)));

mat = zeros(2,i.nstates);
mat(1,intersect(s.I, s.rdt)) = 1;
mat(2,intersect(s.Z, s.rdt)) = 1;
agg.rdt_clin = sparse(mat);

mat = zeros(2,i.nstates);
mat(1,intersect([s.A, s.P, s.I], s.rdt)) = 1;
mat(2,intersect([s.S, s.E, s.R], s.rdt)) = 1;
agg.rdt_comm = sparse(mat);


% --- Counting unnecessary quarantines
mat = zeros(i.nstates);
mat(intersect([s.S, s.E, s.R, s.Z], s.q1), s.q0) = 1;
sel.fpqur = sparse(mat - diag(diag(mat)));

mat = zeros(2, i.nstates);
mat(1, intersect([s.S, s.E, s.R], s.q1)) = 1;
mat(2, intersect(s.Z, s.q1)) = 1;
agg.fpqur = sparse(mat);

% --- Counting all quarantines
mat = zeros(i.nstates);
mat(intersect([s.S, s.E, s.R, s.Z, s.I, s.A, s.P], s.q1), s.q0) = 1;
sel.qur = sparse(mat - diag(diag(mat)));

mat = zeros(2, i.nstates);
mat(1, intersect([s.S, s.E, s.R, s.I, s.A, s.P], s.q1)) = 1;
mat(2, intersect(s.Z, s.q1)) = 1;
agg.qur = sparse(mat);


% The above lines are setting up matrices to help track the 'incidence-like' 
% terms, for example the incidence of symptomatic illness. Don't worry too 
% much about them for now, will explain



% --- Set up paramerers for uncertainty estimates -------------------------
names = {'p_sympto','p_c','p_sympto_chrel'};
lgths = [         1,    1,               1];
xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.nx = lim;

% --- Set up uncertainty bounds -------------------------------------------
bds = zeros(xi.nx,2);
bds(xi.p_sympto,:)       = [0.5, 0.8];
bds(xi.p_c,:)            = [0.5, 1];
bds(xi.p_sympto_chrel,:) = [0.25, 0.75];
% bds(xi.p_NCS,:)      = [0.01, 0.1];
% bds(xi.r_careseek,:) = [0.01, 0.1];
prm.bounds           = bds';



% --- Set up the parameters -----------------------------------------------

R0  = 1;                                                                   % Basic reproduction number
incub_period  = 5;                                                         % Average incubation period (days)
presym_period = 1;                                                         % Average pre-symptomatic period (days)
infec_period  = 5;                                                         % Average infectious period (days)

% --- Set up parameter values
r.incub    = 1/incub_period;                                               % Rate of progression to infectiousness
p.sympto   = 2/3*[0.5 1 1];                                                % Proportion developing symptomatic infection
p.c        = 2/3;                                                          % Relative infectiousness of asymptomatic vs symptomatic infection
p.c2       = p.c;                                                          % Relative infe of pre-symp vs. symp
r.eta      = 1/presym_period;                                              % Rate from pre-symptomatic to symptomatic infection
r.gamma    = 1/infec_period;
% p.hosp     = prop_hosp;
% r.hosp     = 1/5;                                                          % Amongst those with severe disease, average rate of progression to needing hospitalisation
% r.mu       = ??;                                 % Rate of mortality amongst those needing hospitalisation
% r.gamma_h  = (1-mort_on_hosp)*(1/dur_in_hosp);                             % Rate of recovery amongst those needing hospitalisation

% prop_hosp    = [0.04 0.278 0.599];                                         % Age-specific proportions needing hospitalisation (severe disease)
cfr          = [0.01, 0.3, 6.4]/100;                                        % Age-specific case fatality rates
% mort_on_hosp = cfr./prop_hosp;                                             % Amongst those being hospitalised, what proportion would die
% dur_in_hosp  = 10;                                                         % Average days in hospital
r.mu        = cfr./(1-cfr)*r.gamma;

p.imm            = 0;
r.waning         = 0;
prm.PCR_capacity = 1e5;
prm.pNCS         = 0.09;%%%
p.LFA_cs         = 0;
% p.sympto = 0.05;                                                         % What proportion of the population have COVID symptoms, but no COVID

% Interventions
p.sens      = [0.99, 0.8];                                                 % Sensitivity of PCR and RDT testing
p.spec      = [0.99, 0.98];                                                % Specificity of two tests
r.Dx        = [1, 24];                                                     % Turnaround times once test initiated
r.hold      = 2;                                                           % For PCR, default queuing time to initiate test
r.screen    = zeros(3,2);                                                  % Proportion of whole population screened per day
r.careseek  = 0;                                                           % Rate at which symptomatics seek care (for facility-based contact with PCR)
r.Zcareseek = 0;
r.quar      = 1/7;                                                         % Rate of release from quarantine (used only for non-COVID population)

% In the above lines, notice that we use structures once again, i.e. p and
% r, as a coding convenience. In general, all proportions are listed under
% p and all rates are listed under r. Type 'p' and you'll see the values of
% all proportions involved in the model.


% Get India parameters
prm.N       = [7554531	11954901	818671];
prm.connmat = 1;
% prm.mixmat  = [5.59	2.57	0.08
%               2.18	5.56	0.08
%               0.06	0.12	0.01];       
% 
% prm.mixmat  = [0.10	0.486	0.011
%               1.14	6.50	0.14;
%               0.097	0.517	0.016]; 
          
prm.mixmat  = [2.407 0.971 0.052
               1.496 3.681 0.134
               0.093 0.143 0.022]; 
% 'mixmat' is the contact matrix between different age groups. Note: this is
% only provisional, to be updated in light of recent data
           
r.beta = 1;
r.beta = 1/find_R0(p, r, i, s, gps, prm);
% This step is just to ensure that beta is chosen so that R0 = 1, at this 
% preparatory stage. When we come round to simulating the model, we can 
% then simply multiply beta by the desired value of R0.

% find_R0(p, r, i, s, gps, prm)
% r.careseek = 100; prm.r = r;
% find_R0(p, r, i, s, gps, prm)

prm.r = r; prm.p = p;

save Model_setup_Feb22;
