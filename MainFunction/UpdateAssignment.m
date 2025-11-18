function  [T] = UpdateAssignment (T,x,new_t,tmp,prm)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

for k=1:prm.m
    T{k}.Pt_x(x) = new_t;
    % update centroid p(y|t) by merging x to new_t
    T{k}.Py_t = tmp{k}.Py_t;
    mass = tmp{k}.Px+T{k}.Pt(new_t);
    pi1 = tmp{k}.Px/mass;
    pi2 = T{k}.Pt(new_t)/mass;
    T{k}.Py_t (:,new_t) = pi1*tmp{k}.Py_x+pi2*T{k}.Py_t(:,new_t);
    
    % update priors
    T{k}.Pt = tmp{k}.Pt;
    T{k}.Pt(new_t) = mass;
end
end

