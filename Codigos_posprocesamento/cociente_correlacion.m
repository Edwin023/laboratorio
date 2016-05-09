A1=load('correlation4.mat','g1_0','g1_dois','Db','tau');
%A2=load('correlation1.mat','g1_0','g1_dois','Db','tau');
A3=load('correlation3.mat','g1_0','g1_dois','Db','tau');
A4=load('correlation.mat','g1_0','g1_dois','Db','tau');

g1_1=A1.g1_0;
g1_1_duas=A1.g1_dois;

%g1_2=A2.g1_0;
%g1_2_duas=A2.g1_dois;

g1_3=A3.g1_0;
g1_3_duas=A3.g1_dois;

g1_4=A4.g1_0;
g1_4_duas=A4.g1_dois;
figure
semilogx(A1.tau,g1_1(4,:),'bo',A1.tau,g1_1_duas(4,:),'-b',A1.tau,g1_3(4,:),'ro',A1.tau,g1_3_duas(4,:),'-r',A4.tau,g1_4(4,:),'ko',A4.tau,g1_4_duas(4,:),'-k')
axis([1e-7 1e-1 0 1])
legend('Db2/Db1 mc=1','Db2/Db1 Two layer=1','Db2/Db1 mc=100','Db2/Db1 Two layer=100','Db2/Db1 mc=1000','Db2/Db1 Two layer=1000')
 figure
% semilogx(A2.tau,g1_2(4,:),'bo',A2.tau,g1_2_duas(4,:),'-b',A4.tau,g1_4(4,:),'ro',A4.tau,g1_4_duas(4,:),'-r')
% axis([1e-7 1e-1 0 1]) 
% legend('Db2/Db1 mc=10','Db2/Db1 Two layer=10','Db2/Db1 mc=1000','Db2/Db1 Two layer=1000')