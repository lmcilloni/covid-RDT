function R0 = find_R0(p, r, i, s, gps, prm)

M      = make_model_iterative(p, r, i, s, gps, prm);
% M      = make_model_iterative(p, r, i, s, gps, prm);
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
R0 = max(eigs(V\F));