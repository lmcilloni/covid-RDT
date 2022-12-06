function [p,r,prm] = alloc_parameters3(x, p, r, xi, i, s, gps, prm, R0)

% --- Allocate parameters -------------------------------------------------
p.sympto   = x(xi.p_sympto)*[x(xi.p_sympto_chrel) 1 1];
p.c        = x(xi.p_c);
% prm.pNCS   = x(xi.p_NCS);
% r.careseek = x(xi.r_careseek);

% --- Find value of beta required to meet specified value of R0 -----------
r.beta = 1;
M      = make_model_iterative(p, r, i, s, gps, prm);
Iinds  = [s.Eq0, s.infectious];

% Construct F matrix
m = zeros(i.nstates);
m(s.Eq0,:) = repmat(M.lambda,4,1);
% Adjust for initial conditions
vec = [prm.N, zeros(1,9)]';
m(s.Eq0,:) = m(s.Eq0,:).*repmat(vec,1,i.nstates);
F = m(Iinds, Iinds);
% Construct V matrix
V = -M.lin(Iinds, Iinds);
% Bring them together
tmp = max(eigs(V\F));
r.beta = R0/tmp;


% M          = make_model2(p, r, i, s, gps, prm);
% Iinds      = [s.E1, s.infectious];
% % Construct F matrix
% m = zeros(i.nstates);
% m(s.E1,:) = M.lambda;
% % Adjust for initial conditions
% m(s.E1,:) = m(s.E1,:).*repmat(prm.N',1,i.nstates);
% F = m(Iinds, Iinds);
% % Construct V matrix
% V = -M.lin(Iinds, Iinds);
% % Bring them together
% tmp = max(eigs(V\F));
% r.beta = R0/tmp;