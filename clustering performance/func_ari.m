function ARI_value = func_ari(realLabel,predictLabel)
% 生成级联列表
C=Contingency(realLabel,predictLabel);
% addpath 'G:\matlab\Matlab\bin\test-me\higher-order-organization-matlab-master\IB_Release1_0';
% C = confusionMatrix(realLabel,predictLabel);

% 求ARI的各部分数值
n = sum(sum(C));
X = 0.5*sum(sum(C.^2 - C));
ai = sum(C,2);
Y = 0.5*sum(ai.^2 - ai);
bi = sum(C,1);
Z = 0.5*sum(bi.^2 - bi);
exp = Y*Z/nchoosek(n,2);
% 求ARI值
ARI_value = (X-exp)/(0.5*(Y+Z)-exp);

function Cont=Contingency(Mem1,Mem2)
if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1)
       Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1; 
end
end

end

%% 说明：
% 自己写的混淆矩阵的那个代码，实际上就是级联列表的代码，但是对于行列数需要做出修改！  已修改！
