function  [T]= RandomPartition(cellX,cellInp,prm)
%define the structure body T for every output
T=cell(prm.m,1);

for i=1:prm.m
    X = size(cellX{i},1);
    Y = size(cellX{i},2);
    T{i}.size = prm.Tsize;
    T{i}.Pt_x = zeros(1,X);
    T{i}.Pt = zeros(1,prm.Tsize);
    T{i}.Py_t = zeros(Y,prm.Tsize);
    T{i}.costs = zeros(1,X);
    T{i}.Counter = 1;
end

% % Produce a random partition as an initialization point
% the random partition of first view
row=size(cellX{1},1);
perm = randperm(row);
AvgSize =row/prm.Tsize;
i2 = 0;
for t=1:prm.Tsize
i1 = i2+1;
i2 = min(ceil(t*AvgSize),row);
inds = perm(i1:i2);
T{1}.Pt_x(inds) = t;
T{1}.Pt(t) = sum(cellInp{1}.Px(inds));
T{1}.Py_t(:,t) = sum(cellInp{1}.Pyx(:,inds),2)./T{1}.Pt(t);
end

for i=2:prm.m % the rest of the views including shared view
    for t=1:prm.Tsize
        inds =find(T{1}.Pt_x==t);
        T{i}.Pt_x(inds) = t;
        T{i}.Pt(t) = sum(cellInp{i}.Px(inds));
        T{i}.Py_t(:,t) = sum(cellInp{i}.Pyx(:,inds),2)./T{i}.Pt(t);
    end
end

end

