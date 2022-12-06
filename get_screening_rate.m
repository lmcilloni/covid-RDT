% For a given deployment of RDTs in clinics, and for a given overall
% consumption of RDTs, what is the community screening rate needed?
% Function helps us to control for the total consumption of RDTs

function [rscreen, soln3] = get_screening_rate(p, r, nrdt, init, tf, s, i, agg, sel, gps, prm)

% Now, find the screening rate needed to give the same consumption of RDTs
hi = 0; lo = 1; nrdt2 = 10*nrdt;

while abs(1-nrdt2/nrdt)>1e-6
    rscreen = mean([lo, hi]);
    
    % PCR testing in clinical, plus community screening with RDT
    p3 = p; r3 = r;
    r3.screen = [0 rscreen];
    M3 = make_model_iterative(p3, r3, i, s, gps, prm);
    
    geq = @(t,in) goveqs_basis(t, in, M3, i, s, r3, prm, agg, sel);
    [t,soln3] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));

    nrdt2 = sum(soln3(end,i.aux.rdt));
    if nrdt2>nrdt
        lo = rscreen;
    else
        hi = rscreen;
    end
end

