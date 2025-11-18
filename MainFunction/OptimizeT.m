function  [T] = OptimizeT (T,cellX,cellInp,prm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

LoopCounter = 1;
NoChange =0;
done=0;

%% intial weights
Weights=ones(prm.m,1)*(1/(prm.m));

while ~done
    change=0;
    %% updating T
    for x=1:size(cellX{1},1)                % documents in every dataset should be same
        [tmp,t,inds_t]=ProcessTmp(cellInp,T,prm,x);
        if length(inds_t)>1
            [tmp] =Reduce_x (tmp,t,inds_t,x,cellInp,prm);  %change every Pt and P(y,t)
            [cellCosts] =MergeCosts (tmp,prm);
            [~,new_t] =MinSumCosts(cellCosts,prm,Weights);
            if new_t==t
                NoChange= NoChange+1;
            else
                NoChange= 0;
                change= change+1;
                [T] =UpdateAssignment (T,x,new_t,tmp,prm);
            end
        end  %end if length
    end %end for x=1:
    [T, done] =CheckConvergence (T,prm,LoopCounter,change);
    LoopCounter= LoopCounter+1;
    
    %% Update weights W:
    % computing objective
    obj=zeros(prm.m,1);
    ITY=zeros(prm.m,1);
    HT=zeros(prm.m,1);
    for i=1:prm.m
        dd=T{i}.Pt ;
        PPt = repmat(dd,size(cellX{i},2),1);   % PPT is used to calculate P(t,y) in the next step
        Pty = T{i}.Py_t.*PPt;
        [aa,~,cc] = MTC_local_MI (Pty,dd);
        ITY(i)=aa;
        HT(i)=cc;
        % 此处其实有个Q：在context和content中为什么HT是m/2,只取了前一半的HT？

    end
    
    temp_denominator=0;
    for i=1:prm.m
        obj(i)=HT(i)-(1/prm.inv_beta)*ITY(i);
        temp_denominator = temp_denominator+exp(-obj(i)/prm.theta);
    end
    
%  这一部分是权重值的输出,可能需要考虑一下如何输出最终的值,实在不行,全部输出,自己找最后一个 %%
    for i=1:prm.m
        if temp_denominator==0
            Weights(i)=0.5;
        else
            Weights(i) = exp(-obj(i)/prm.theta)/temp_denominator;  % AMCIB  AMCIBW的话修改此处的权值
              % disp([' Weights' num2str(i) ' is:  ' num2str(Weights(i))]);
%             % Debug*  %% 输出的时候，可以参照上面的这行代码 %%
        end 
    end
    
end  %end while


% this also needed to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for i=1:prm.m
    dd=T{i}.Pt ;
    PPt = repmat(dd,size(cellX{i},2),1);   % PPT is used to calculate P(t,y) in the next step
    Pty = T{i}.Py_t.*PPt;
    [aa,bb,cc] = MTC_local_MI (Pty,dd);
    T{i}.Ity=aa;
    T{i}.Hy=bb;
    T{i}.Ht=cc;
    T{i}.L = T{i}.Ity;
%     T{i}.L = T{i}.Ity - prm.inv_beta*T{i}.Ht;
    % There seems to be a saying that I (T;X) = Ht
end

end
