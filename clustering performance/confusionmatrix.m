function conM = confusionmatrix(indx,L,clusters,categories)
%clusters 是聚类后的种类数，category是原来几个类别
%index是聚类结果，L是data.L1类标签

conM = zeros(categories,clusters);
u = unique(indx);  
%取集合a的不重复 元素 构成的向量
for i = 1:length(indx)  
     %length(x)表示向量x的坐标分量的个数，即向量的维数
     %返回行数和列数两者之间数值较大的那一个值，即相当于执行了max(size(A))
    r = find(L(i,:) ~= 0);
    %MATLAB中"："的使用非常灵活，简单的说就是代表矩阵中的全行或全列元素
    %冒号表示所有的index，因此Y(i,j)表示Y的第i行第j列，Y(:,i)就表示整个第i列了
    c = find(u == indx(i));
    %index中第i 个元素
    %满足find中括号里面的条件的元素位置
    conM(r,c) = conM(r,c) + 1;
end