%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MI.m
%
% IBA_1.0 - Information Bottleneck Algorithms for "hard" clustering, Version 1.0
% Copyright (C) 2003 Noam Slonim, Nir Friedman, and Naftali Tishby
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I,Hx,Hy] = MI(P)
% MI(P) returns the mutual information and the entropies of the joint distribution P, where:
%   I(p(x,y)) = sum_{x,y} p(x,y)*log(p(x,y)/p(x)p(y));

if min(size(P))<2
  warning ('This is not a JOINT distribution');
  Hx = Ent(P);  % Ent returns the Entropy of p, given by -p*log(p) where p is a normalized distribution
  Hy = 0;
  I = 0;
end

if min(min(P)) < 0,
  error ('Negative probability values');
end;

if (abs(sum(sum(P))-1) >= 1e-6)
  fprintf ('sum_{x,y}P(x,y)=%f\n',full(sum(sum(P))));
  error ('Unnormalized joint distribution');
end;

Px=sum(P');
Hx=Ent(Px);
Py=sum(P);
Hy=Ent(Py);

inds = find(P);
I = Hx + Hy + sum(sum(P(inds).*MyLog(P(inds))));	

return
