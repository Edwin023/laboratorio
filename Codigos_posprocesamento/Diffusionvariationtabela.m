function Diffusionvariationtabela(num)
% clc
% close all
% %clear
load('sample.mat');
load('tau.mat');
%Db=[1e-8 1e-5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Db2=[1.0e-6 1.5e-6 3e-6];
s.g1=zeros(length(Db2),length(taufit));

for i=1:length(Db2)

Db=[1.5e-8 Db2(i)];
%tau=1e-7:1e-6:1e-0;
lambda=852e-7;

[~,g1_mc,~]=correlacion_p_mn('sample',lambda,Db,taufit,1,1,0);

s(i).g1=g1_mc;
clear g1;
end
save(['phantom_' num2str(num) 'mm_b1.mat'])
clear s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Db2=[1.0e-6 1e-6 1.5e-6 3e-6];
Db1=[1e-8 1.5e-8 1.5e-8 1.5e-8];

for i=1:length(Db2)

Db=[Db1(i) Db2(i)];
%tau=1e-7:1e-6:1e-0;
lambda=852e-7;

[~,g1_mc,~]=correlacion_p_mn('sample',lambda,Db,taufit,1,1,0);

s(i).g1=g1_mc;
clear g1;
end
save(['phantom_' num2str(num) 'mm_b2.mat'])
clear s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
Db2=[1.0e-6 1.5e-6 3e-6];
%s.g1=zeros(length(Db2),length(taufit));
Db1=[1e-8 1.5e-8 3e-8];

for i=1:length(Db2)

Db=[Db1(i) Db2(i)];
%tau=1e-7:1e-6:1e-0;
lambda=852e-7;

[~,g1_mc,~]=correlacion_p_mn('sample',lambda,Db,taufit,1,1,0);

s(i).g1=g1_mc;
clear g1;
end
save(['phantom_' num2str(num) 'mm_b3.mat'])
clear s;



%%Cuarta simulacion
Db2=[1.0e-6 1.25e-6 1.5e-6 2e-6];
Db1=[1e-8 1.1e-8 1.25e-8 1.5e-8];
for i=1:length(Db2)

Db=[Db1(i) Db2(i)];
%tau=1e-7:1e-6:1e-0;
lambda=852e-7;

[~,g1_mc,~]=correlacion_p_mn('sample',lambda,Db,taufit,1,1,0);

s(i).g1=g1_mc;
clear g1;
end
save(['phantom_' num2str(num) 'mm_b4.mat'])
clear s;



%%intensidas
[~,I]=mcxyz_read_his('sample',1);
fid=fopen('sample1_F.bin','rb');
fluencia=fread(fid,'float32');
fluencia=reshape(fluencia,[240 240 153]);
fclose(fid);

%%Read Total Photon simulated
fid=fopen('photon_number.txt','r');
NTphotons=fscanf(fid,'%f',1);
fclose(fid);

save(['phantom_' num2str(num) 'mm_e-10.mat'])
clear fluencia;
end
