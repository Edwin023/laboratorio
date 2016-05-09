clc
close all
clear
load('sample.mat');
load('tau.mat');
%Db=[1e-8 1e-5];

Db2=[1.5e-8 2.0e-8 2.5e-8 3.0e-8 3.5e-8];
s.g1=zeros(length(Db2),length(taufit));
for i=1:length(Db2)

Db=[1e-10 Db2(i)];
%tau=1e-7:1e-6:1e-0;
lambda=852e-7;

[nPhotons2,g1_mc,mc2]=correlacion_p_mn('sample',lambda,Db,taufit,1,1,0);

s(i).g1=g1_mc;
clear g1;
end
[~,I]=mcxyz_read_his('sample',1);
save('phantom_8mm_e-10')