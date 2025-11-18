%% comparative matrix
% Idx=(CR(1,:))'; %sIB中，每次取得是迭代次数的每一行，其代表每次迭代的结果
Idx=CR;%Idx是最终聚类结果，维度为1*N
D=gnd;%gnd是真实的数据，维度为N*1

%% AC evaluation      
 %CL = bestMap(D,Idx');
   AC = averageAC(Idx,D) ;
   
 %% MI evaluation
 num=size(CR,1);
 MIhat=zeros(1,num);
 for i=1:num
    MIhat(i)= MutualInfo(D,CR(i,:)); % pay attention to the CR(i,:);
 end
 MIhat=sum(MIhat,2)/num;
    
        
 %% 生成混淆矩阵 
 %conM = confusionmatrix(indx,L,clusters,categories)
 %clusters 是聚类后的种类数，category是原来几个类别
 %index是聚类结果，L是data.L1类标签
 conM = confusionmatrix(Idx',D,2,2);
           
 %% 计算精确度评估结果
 %Pre = sum(max(conM))/length(CL);
        