
function [kl] = KL(p,q,Chk_q)
% KL returns the Kullbac-Leibler divergence, where:
%   KL(p,q) = p*log(p/q) 

inds = find(p);
p = p(inds);
q = q(inds);

if Chk_q % not necessary when called from JS
    if min(q)==0
      kl=inf;
      return
    end
end

kl = sum(p.*MyLog(p./q));

return
