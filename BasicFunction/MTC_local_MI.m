function [I,Hy,Ht] = MTC_local_MI (P,Pt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

Py=sum(P');
Hy = -sum(Py.*MyLog(Py));
Ht = -sum(Pt.*MyLog(Pt));
inds = find(P);
I = Hy + Ht + sum(sum(P(inds).*MyLog(P(inds))) );	% I = Hy+Ht-Hyt


end

