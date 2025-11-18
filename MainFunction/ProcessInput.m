function  [cellInp]=ProcessInput(cellX,prm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
cellInp=cell(prm.m,1);

 for i=1:prm.m
 if ((min(sum(cellX{i},1))<=0) | (min(sum(cellX{i},2))<=0))
  error ('zero priors in M - empty columns or rows are not allowed')
 end
 
 if prm.UniformPrior
   cellX{i}=NormMat(cellX{i});
 end
 
cellInp{i}.Py_x = (NormMat(cellX{i}))';
cellInp{i}.Pyx = cellX{i}'./sum(sum(cellX{i}));
cellInp{i}.Px = sum(cellInp{i}.Pyx,1);
cellInp{i}.Py = sum(cellInp{i}.Pyx,2);

[cellInp{i}.Ixy] = MI(cellInp{i}.Pyx);
 end 

end

