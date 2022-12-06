function M = make_model_iterative_dynamic(p, r, i, s, gps, prm)

% Every model can be fully specified by two matrices: one capturing all
% linear transitions, and the other capturing the (nonlinear) infection
% process. Writing models in this way can be very helpful for models with
% large numbers of states, expecially when a basic model is stratified
% (e.g. in our case, a basic structure being stratified by age).
% 
% -------------------------------------------------------------------------
% --- First, set up the matrix to represent all linear (i.e. non-infection-
% related) transitions

m  = zeros(i.nstates);
m1 = zeros(i.nstates);                                                     % Separate matrix for rate of diagnosis with PCR
m2 = zeros(i.nstates);
for ia = 1:length(gps.age)
    age = gps.age{ia};
    
    for iq = 1:length(gps.quar)
        quar = gps.quar{iq};
        
        for it = 1:length(gps.test)
            test = gps.test{it};
            
            getind = @(st) i.(st).(test).(quar).(age);
            S = getind('S');
            E = getind('E');
            A = getind('A');
            P = getind('P');
            I = getind('I');
            R = getind('R');
%             H = i.H.(age);
            D = i.D.(age);
            
%             Incubation
            source  = E;
            destins =                [A,                P];
            rates   = [(1-p.sympto(ia)),     p.sympto(ia)]*r.incub;
            m(destins, source) = m(destins, source) + rates';
            
%             Pre-symptomatic to symptomatic
            source = P;
            destin = I;
            rates  = r.eta;
            m(destin, source) = m(destin, source) + rates';
            
%             % Hospitalization
%             source = I_stub.svr;
%             destin = H;
%             rate   = r.hosp;
%             m(destin, source) = m(destin, source) + rate;
            
%             Recovery
            sources = [A, I];
            destin  = R;
            rate    = r.gamma;
            m(destin, sources) = m(destin, sources) + rate;
            
        end
    end
    
%     --- Mortality
    source  = I;
    destin  = D;
    rate    = r.mu(ia);
    m(destin, source) = m(destin, source) + rate;
    
%     % --- Recovery amongst hospitalised
%     source  = H;
%     destin  = i.R.null.q0.(age);
%     rate    = r.gamma_h(ia);
%     m(destin, source) = m(destin, source) + rate;
    
end

% --- Careseeking amongst those with symptomatic COVID
sources = intersect(intersect(s.I, s.q0),s.null);
destins = intersect(intersect(s.I, s.q0),s.pcra);
rates   = r.careseek*(1-p.LFA_cs);
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

sources  = intersect(intersect(s.I, s.q0),s.null);
destins  = intersect(intersect(s.I, s.q0),s.rdt);
rates    = r.careseek*p.LFA_cs;
inds     = sub2ind([i.nstates,i.nstates], destins, sources);
m2(inds) = m2(inds) + rates;

sources  = intersect(intersect(s.I, s.q0),s.rdt);
destins  = intersect(intersect(s.I, s.q0),s.pcra);
rates    = r.Dx(2)*p.sens(2); % pcr confirmation (true pos)
inds     = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

sources  = intersect(intersect(s.I, s.q0),s.rdt);
destins  = intersect(intersect(s.I, s.q0),s.null);
rates    = r.Dx(2)*(1-p.sens(2)); % pcr conf (false neg)
inds     = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;


% DETECTING true positives (of PCR conf)
sources = intersect(intersect(s.pcrb,s.q0), s.I);
destins = intersect(intersect(s.null,s.q1), s.I);
rates   = r.Dx(1)*p.sens(1);
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

% MISSING true positives (of PCR conf)
sources = intersect(intersect(s.pcrb,s.q0), s.I);
destins = intersect(intersect(s.null,s.q0), s.I);
rates   = r.Dx(1)*(1-p.sens(1));
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;


% --- Careseeking amongst those with symptoms but no COVID
sources = intersect(intersect(s.Z, s.q0),s.null);
destins = intersect(intersect(s.Z, s.q0),s.pcra);
rates   = r.Zcareseek*(1-p.LFA_cs);
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

sources  = intersect(intersect(s.Z, s.q0),s.null);
destins  = intersect(intersect(s.Z, s.q0),s.rdt);
rates    = r.Zcareseek*p.LFA_cs;
inds     = sub2ind([i.nstates,i.nstates], destins, sources);
m2(inds) = m2(inds) + rates;

sources = intersect(intersect(s.Z, s.q0),s.rdt);
destins = intersect(intersect(s.Z, s.q0),s.pcra);
rates   = r.Dx(2)*(1-p.spec(2)); % pcr confirmation (false pos)
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

sources = intersect(intersect(s.Z, s.q0),s.rdt);
destins = intersect(intersect(s.Z, s.q0),s.null);
rates   = r.Dx(2)*p.spec(2); % pcr confirmation (true neg)
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

% MISDIAGNOSING false positives (of PCR conf)
sources = intersect(intersect(s.pcrb,s.q0),s.Z);
destins = intersect(intersect(s.null,s.q1),s.Z);
rates   = r.Dx(1)*(1-p.spec(1));
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

% IDENTIFYING true negatives (of PCR conf)
sources = intersect(intersect(s.pcrb,s.q0),s.Z);
destins = intersect(intersect(s.null,s.q0),s.Z);
rates   = r.Dx(1)*p.spec(1);
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;


% --- Release from quarantine
sources = s.q1;
destins = s.q0;
rates   = r.quar;
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m(inds) = m(inds) + rates;

% --- Population screening and test results -------------------------------

% test1s = {'pcra','rdt'}; test2s = {'pcrb','rdt'};
% 
% for id = 1:2
%     st1 = s.(test1s{id});
%     st2 = s.(test2s{id});
    
    for ia = 1:length(gps.age)
        age = gps.age{ia};        
%         Initiating testing through screening
        sources = intersect(intersect(intersect(s.null,s.q0),s.main),s.(age)); % NB: screening does not apply to Z, since non-Covid compartments in main population provide basis for false-positives
        destins = intersect(intersect(intersect(s.rdt,s.q0),s.main),s.(age));
        rates   = r.screen(ia,2);
        inds    = sub2ind([i.nstates,i.nstates], destins, sources);
        m(inds) = m(inds) + rates;          
    end
    
%     DETECTING true positives
    sources = intersect(intersect(s.rdt,s.q0),[s.I, s.A, s.P]);
    destins = intersect(intersect(s.pcra,s.q0),[s.I, s.A, s.P]);
    rates   = r.Dx(2)*p.sens(2);
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     MISSING true positives
    sources = intersect(intersect(s.rdt,s.q0),[s.I, s.A, s.P]);
    destins = intersect(intersect(s.null,s.q0),[s.I, s.A, s.P]);
    rates   = r.Dx(2)*(1-p.sens(2));
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     MISDIAGNOSING false positives
    sources = intersect(intersect(s.rdt,s.q0),[s.S, s.E, s.R]);
    destins = intersect(intersect(s.pcra,s.q0),[s.S, s.E, s.R]);
    rates   = r.Dx(2)*(1-p.spec(2));
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     IDENTIFYING true negatives
    sources = intersect(intersect(s.rdt,s.q0),[s.S, s.E, s.R]);
    destins = intersect(intersect(s.null,s.q0),[s.S, s.E, s.R]);
    rates   = r.Dx(2)*p.spec(2);
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;

    
    %%% test 2
%     DETECTING true positives
    sources = intersect(intersect(s.pcrb,s.q0),[s.I, s.A, s.P]);
    destins = intersect(intersect(s.null,s.q1),[s.I, s.A, s.P]); % isolation
    rates   = r.Dx(1)*p.sens(1);
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     MISSING true positives
    sources = intersect(intersect(s.pcrb,s.q0),[s.I, s.A, s.P]);
    destins = intersect(intersect(s.null,s.q0),[s.I, s.A, s.P]);
    rates   = r.Dx(1)*(1-p.sens(1));
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     MISDIAGNOSING false positives
    sources = intersect(intersect(s.pcrb,s.q0),[s.S, s.E, s.R]);
    destins = intersect(intersect(s.null,s.q1),[s.S, s.E, s.R]); % isolation
    rates   = r.Dx(1)*(1-p.spec(1));
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
%     IDENTIFYING true negatives
    sources = intersect(intersect(s.pcrb,s.q0),[s.S, s.E, s.R]);
    destins = intersect(intersect(s.null,s.q0),[s.S, s.E, s.R]);
    rates   = r.Dx(1)*p.spec(1);
    inds    = sub2ind([i.nstates,i.nstates], destins, sources);
    m(inds) = m(inds) + rates;
    
% end

% --- PCR testing, delay due to carrying capacity (separate matrix)
sources = intersect(s.pcra,s.q0);
destins = intersect(s.pcrb,s.q0);
rates   = r.hold;
inds    = sub2ind([i.nstates,i.nstates], destins, sources);
m1(inds) = m1(inds) + rates;


% --- Bring them all together ---------------------------------------------
M.lin           = sparse(m - diag(sum(m,1)));
M.lin_PCR_queue = sparse(m1 - diag(sum(m1,1)));
M.lin_RDT_clin  = sparse(m2 - diag(sum(m2,1)));

% --- Nonlinear (infection-related) matrix --------------------------------
% 
% Finally, set up the infection matrix, used to capture the
% nonlinear transition process. When this matrix multiplies the state
% vector, it gives the force-of-infection for each of the 3 age groups -
% hence below, it's labelled 'lambda'
% 
% keyboard;
% 
% Adjust contact matrix by population numbers to get a repeating block
rpt = prm.mixmat./repmat(prm.N,length(gps.age),1);
% Allocate these terms to all infectious compartments
m = zeros(length(gps.age),i.nstates);

cols = intersect(s.A,s.infectious);
m(:,cols) = repmat(rpt,1,length(cols)/length(gps.age))*p.c;                % Also adjusting for infectivity of a-pre-symptomatics

cols = intersect(s.P,s.infectious);
m(:,cols) = repmat(rpt,1,length(cols)/length(gps.age))*p.c2;               % for sens analysis we change pre-symp infectivity

cols = intersect(s.I,s.infectious);
m(:,cols) = repmat(rpt,1,length(cols)/length(gps.age));


% cols = intersect(s.infectious,s.svr);
% m(:,cols) = repmat(rpt,1,length(cols)/length(gps.age));

M.lambda = sparse(m)*r.beta;

