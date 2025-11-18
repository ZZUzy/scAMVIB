
function [H] = Ent(p)
% Ent returns the Entropy of p, given by -p*log(p) where p is a normalized distribution

if min(p) < 0,
  error ('Negative probability values');
end;

inds = find(p);
p = p(inds);

if (abs(sum(p)-1) >= 1e-6),
    err_str = sprintf ('Unnormalized distribution: sum_{x}P(x)=%f\n',full(sum(p)));
    error (err_str);
end;

H = -sum(p.*MyLog(p));

return;

	
