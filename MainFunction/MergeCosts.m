function [cellCosts] =MergeCosts (tmp,prm)
% Costs(t) = (Px+Pt(t))*JS(Py_x,Py_t(:,t)) - inv_beta*H([pi1 pi2])

cellCosts = cell(prm.m,1);
for i=1:prm.m
    cellCosts{i}=zeros(1,prm.Tsize);
end

for i=1:prm.m
        for t=1:prm.Tsize
            Pnew= tmp{i}.Px+tmp{i}.Pt(t);
            pi1= tmp{i}.Px/Pnew;
            pi2= tmp{i}.Pt(t)/Pnew;
            tmpJS= JS(tmp{i}.Py_x,tmp{i}.Py_t(:,t),pi1,pi2);
            if prm.inv_beta~=0
                cellCosts{i}(t) = Pnew*(tmpJS - prm.inv_beta*Ent([pi1 pi2]));
            else
                cellCosts{i}(t) = Pnew*(tmpJS);
            end
        end
end
 
end
