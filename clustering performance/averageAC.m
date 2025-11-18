function AC = averageAC(CR,gnd)  

%% 对多种聚类结果求平均值可以减小误差,同时也可用方差发现这几个平均值间的波动大小。
%CR表示聚类结果，gnd表示真实的源数据;
T = size(CR,1);
AC = zeros(1,T);

for i = 1:T
    res = bestMap(gnd,CR(i,:));
    AC(i) = length(find(gnd == res))/length(gnd);
end

%% get the average
%  AC=sum(AC,2)/T;


    