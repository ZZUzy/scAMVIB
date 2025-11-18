function Results = evaluation(rpath)

%files = char('Binary_1.mat', 'Binary_2.mat' ,'Binary_3.mat', 'Multi5_1.mat','Multi5_2.mat' ,'Multi5_3.mat', 'Multi10_1.mat', 'Multi10_2.mat' ,'Multi10_3.mat');

files = char('Multi5_1.mat','Multi5_2.mat' ,'Multi5_3.mat', 'Multi10_1.mat', 'Multi10_2.mat' ,'Multi10_3.mat');
 %多个数据集产生的聚类结果的文件集合，所有的idx放一块等于一块循环处理了，不用一个一个处理了！

for i = 1:size(files,1)
    name = deblank(files(i,:));
    in_path = fullfile(rpath,name);
    load(in_path);
    s = ['L =',name(1:end-4),'.L;'];
    eval(s);
    s = ['sMIBCC =',name(1:end-4),'.sMIBCC;'];
    eval(s);
    s = ['sIB =',name(1:end-4),'.sIB;'];
    eval(s)
%     s = ['ITCC =',name(1:end-4),'.ITCC;'];
%     eval(s)
    s = ['SyIB =',name(1:end-4),'.SyIB;'];
    eval(s)
    clear files(i).name(1:end-4)
    [value TL] = max(L,[],2);
    AC = zeros(4,12);
    Pre = zeros(4,12);
    NMI = zeros(4,12);
    for j = 1:10
        CL = bestMap(TL,sIB.Ptx_x(j,:));
        AC(1,j) = length(find(TL == CL))/length(CL); 
        NMI(1,j) = MutualInfo(TL,CL);
        conM = confusionmatrix(sIB.Ptx_x(j,:),L,size(L,2),size(L,2));
        Pre(1,j) = sum(max(conM))/length(CL);
        
        CL = bestMap(TL,sMIBCC.Txty(j).Pt_x);
        AC(2,j) = length(find(TL == CL))/length(CL);
        NMI(2,j) = MutualInfo(TL,CL);
        conM = confusionmatrix(sMIBCC.Txty(j).Pt_x,L,size(L,2),size(L,2));
        Pre(2,j) = sum(max(conM))/length(CL);
        
%         CL = bestMap(TL,ITCC(j,:));
%         AC(3,j) = length(find(TL == CL))/length(CL);
%         NMI(3,j) = MutualInfo(TL,CL);
%         conM = confusionmatrix(ITCC(j,:),L,size(L,2),size(L,2));
%         Pre(3,j) = sum(max(conM))/length(CL);
        
        CL = bestMap(TL,SyIB.Txty(j).Pt_x);
        AC(4,j) = length(find(TL == CL))/length(CL);
        NMI(4,j) = MutualInfo(TL,CL); 
        conM = confusionmatrix(SyIB.Txty(j).Pt_x,L,size(L,2),size(L,2));
        Pre(4,j) = sum(max(conM))/length(CL);
    end
    Pre(1,11) = sum(Pre(1,1:10))/10;
    Pre(1,12) = std(Pre(1,1:10));
    Pre(2,11) = sum(Pre(2,1:10))/10;
    Pre(2,12) = std(Pre(2,1:10));
    Pre(3,11) = sum(Pre(3,1:10))/10;
    Pre(3,12) = std(Pre(3,1:10));
    Pre(4,11) = sum(Pre(4,1:10))/10;
    Pre(4,12) = std(Pre(4,1:10));
    
    AC(1,11) = sum(AC(1,1:10))/10;
    AC(1,12) = std(AC(1,1:10));
    AC(2,11) = sum(AC(2,1:10))/10;
    AC(2,12) = std(AC(2,1:10));
    AC(3,11) = sum(AC(3,1:10))/10;
    AC(3,12) = std(AC(3,1:10));
    AC(4,11) = sum(AC(4,1:10))/10;
    AC(4,12) = std(AC(4,1:10));
    
    NMI(1,11) = sum(NMI(1,1:10))/10;
    NMI(1,12) = std(NMI(1,1:10));
    NMI(2,11) = sum(NMI(2,1:10))/10;
    NMI(2,12) = std(NMI(2,1:10));
    NMI(3,11) = sum(NMI(3,1:10))/10;
    NMI(3,12) = std(NMI(3,1:10));
    NMI(4,11) = sum(NMI(4,1:10))/10;
    NMI(4,12) = std(NMI(4,1:10));
    
    s = ['Results.',name(1:end-4),'.Pre = Pre;'];
    eval(s);
    s = ['Results.',name(1:end-4),'.AC = AC;'];
    eval(s);
    s = ['Results.',name(1:end-4),'.NMI = NMI;'];
    eval(s);
end