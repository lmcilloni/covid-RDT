% For a given deployment of RDTs in clinics, and for a given overall
% consumption of RDTs, what is the community screening rate needed?
% Function helps us to control for the total consumption of RDTs

function [relative_rdt, nrdt2, soln3] = get_num_rdts(rscreen, p, r, nrdt, init, tf, s, i, agg, sel, gps, prm)

% PCR testing in clinical, plus community screening with RDT
p3 = p; r3 = r;
r3.screen = [zeros(3,1), rscreen*ones(3,1)];
M3 = make_model_iterative(p3, r3, i, s, gps, prm);

geq = @(t,in) goveqs_basis(t, in, M3, i, s, r3, prm, agg, sel);
[t,soln3] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));

nrdt2 = sum(soln3(end,i.aux.rdt));
relative_rdt = (1 - nrdt2/nrdt)^2;