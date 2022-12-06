function out = goveqs_basis(t, in, M, i, s, r, prm, agg, sel)

invec = in(1:i.nstates);
out   = zeros(i.nx,1);

% --- Conctruct the linear terms
all_lin = M.lin + M.lin_RDT_clin + (1 - (sum(invec(s.pcrb)))/prm.PCR_capacity)*M.lin_PCR_queue;
out(1:i.nstates) = all_lin*invec;

% Get new infections
lam = M.lambda*invec;
newinfs = repmat(lam,4,1).*invec(s.Sq0);

out(s.Sq0) = out(s.Sq0) - newinfs(:);
out(s.Eq0) = out(s.Eq0) + newinfs(:);

% Get the auxiliaries - these are the 'incidence-like' quantities mentioned in Setup_model being 'counted', e.g. the
% symptomatic incidence. 
out(i.aux.inc)      = agg.inc*(sel.inc.*all_lin)*invec;
% out(i.aux.hosp)     = agg.hosp*(sel.hosp.*all_lin)*invec;
out(i.aux.pcr)      = agg.pcr*(sel.pcr.*all_lin)*invec;
out(i.aux.rdt_comm) = agg.rdt_comm*(sel.rdt.*M.lin)*invec;
out(i.aux.rdt_clin) = agg.rdt_clin*(sel.rdt.*M.lin_RDT_clin)*invec;
out(i.aux.fpqur)    = agg.fpqur*(sel.fpqur.*all_lin)*invec;
out(i.aux.qur)      = agg.qur*(sel.qur.*all_lin)*invec;