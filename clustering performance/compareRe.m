function CPR = compareRe(R)
% dataname = char('Binary_1', 'Binary_2' ,'Binary_3', 'Multi5_1','Multi5_2' ,'Multi5_3', 'Multi10_1', 'Multi10_2' ,'Multi10_3');

dataname = char('Multi5_1','Multi5_2' ,'Multi5_3', 'Multi10_1', 'Multi10_2' ,'Multi10_3');


CPR.AC = zeros(size(dataname,1),4);
CPR.AC_std = zeros(size(dataname,1),4);
CPR.Pre = zeros(size(dataname,1),4);
CPR.Pre_std = zeros(size(dataname,1),4);
CPR.NMI = zeros(size(dataname,1),4);
CPR.NMI_std = zeros(size(dataname,1),4);
for i = 1:size(dataname,1)
    s = ['CPR.AC(',int2str(i),',1)=R.',dataname(i,:),'.AC(1,11);'];
    eval(s);
    s = ['CPR.AC(',int2str(i),',2)=R.',dataname(i,:),'.AC(2,11);'];
    eval(s);
    s = ['CPR.AC(',int2str(i),',3)=R.',dataname(i,:),'.AC(3,11);'];
    eval(s);
    s = ['CPR.AC(',int2str(i),',4)=R.',dataname(i,:),'.AC(4,11);'];
    eval(s);
    
    s = ['CPR.Pre(',int2str(i),',1)=R.',dataname(i,:),'.Pre(1,11);'];
    eval(s);
    s = ['CPR.Pre(',int2str(i),',2)=R.',dataname(i,:),'.Pre(2,11);'];
    eval(s);
    s = ['CPR.Pre(',int2str(i),',3)=R.',dataname(i,:),'.Pre(3,11);'];
    eval(s);
    s = ['CPR.Pre(',int2str(i),',4)=R.',dataname(i,:),'.Pre(4,11);'];
    eval(s);
    
    s = ['CPR.NMI(',int2str(i),',1)=R.',dataname(i,:),'.NMI(1,11);'];
    eval(s);
    s = ['CPR.NMI(',int2str(i),',2)=R.',dataname(i,:),'.NMI(2,11);'];
    eval(s);
    s = ['CPR.NMI(',int2str(i),',3)=R.',dataname(i,:),'.NMI(3,11);'];
    eval(s);
    s = ['CPR.NMI(',int2str(i),',4)=R.',dataname(i,:),'.NMI(4,11);'];
    eval(s);
    
    
    s = ['CPR.AC_std(',int2str(i),',1)=R.',dataname(i,:),'.AC(1,12);'];
    eval(s);
    s = ['CPR.AC_std(',int2str(i),',2)=R.',dataname(i,:),'.AC(2,12);'];
    eval(s);
    s = ['CPR.AC_std(',int2str(i),',3)=R.',dataname(i,:),'.AC(3,12);'];
    eval(s);
    s = ['CPR.AC_std(',int2str(i),',4)=R.',dataname(i,:),'.AC(4,12);'];
    eval(s);
    
    s = ['CPR.Pre_std(',int2str(i),',1)=R.',dataname(i,:),'.Pre(1,12);'];
    eval(s);
    s = ['CPR.Pre_std(',int2str(i),',2)=R.',dataname(i,:),'.Pre(2,12);'];
    eval(s);
    s = ['CPR.Pre_std(',int2str(i),',3)=R.',dataname(i,:),'.Pre(3,12);'];
    eval(s);
    s = ['CPR.Pre_std(',int2str(i),',4)=R.',dataname(i,:),'.Pre(4,12);'];
    eval(s);
    
    s = ['CPR.NMI_std(',int2str(i),',1)=R.',dataname(i,:),'.NMI(1,12);'];
    eval(s);
    s = ['CPR.NMI_std(',int2str(i),',2)=R.',dataname(i,:),'.NMI(2,12);'];
    eval(s);
    s = ['CPR.NMI_std(',int2str(i),',3)=R.',dataname(i,:),'.NMI(3,12);'];
    eval(s);
    s = ['CPR.NMI_std(',int2str(i),',4)=R.',dataname(i,:),'.NMI(4,12);'];
    eval(s);
end