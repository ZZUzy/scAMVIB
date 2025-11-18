
function [js] = JS(p,q,pi1,pi2)
% JS returns the Jensen-Shannon divergence, where:
%   JS(p,q) = pi1*KL(p,ave) + pi2*KL(q,ave), where pi1,pi2>=0, pi1+pi2=1, ave = pi1*p+pi2*q;

if min(pi1,pi2)<=0     
    warning ('zero or negative weights in JS calculation');
    js=0;
    return;
end

ave=pi1*p+pi2*q;

kl1 = KL(p,ave,0);
kl2 = KL(q,ave,0);

js = pi1*kl1 + pi2*kl2;

return
