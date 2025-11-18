function  [tmp,t,inds_t]=ProcessTmp(cellInp,T,prm,x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% process every tmp before optimizing
tmp=cell(prm.m,1);

for i=1:prm.m
    tmp{i}.Py_x = cellInp{i}.Py_x(:,x);
    tmp{i}.Px = cellInp{i}.Px(x);
    tmp{i}.Py_t = T{i}.Py_t;
    tmp{i}.Pt = T{i}.Pt;
    
    t= T{i}.Pt_x(x);
    inds_t = find(T{i}.Pt_x==t);
end

end

