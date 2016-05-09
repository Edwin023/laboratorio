
close all
clear all
clc
% 
% 
% 
% % tiss_prop(1).scattering = 1.010;
% % tiss_prop(1).absorption = 0.0100;
% % tiss_prop(1).refraction = 1.35;
% % tiss_prop(1).anisotropy = 0.001;
% % optpos = [50 50 10;  50 70 10;  50 30 10;  70 50 10;  30 50 10];
% % tMCimg_ejecutor('sample', 100000, [0 5.0e-9 5.0e-9], [100 100 100], optpos, tiss_prop, 1,'windows');
% 
 filename='mc_semi_hist';
  mkdir (filename);
  copyfile('*.m',filename);
  copyfile('*.mat',filename);
  copyfile('mcxyz_mod1.exe',filename);
  cd (filename)

% 
% 
tiss_prop.musv = 10;
tiss_prop.muav = 0.05;
tiss_prop.gv = 0;
% tiss_prop(1,2).musv = 10;
% tiss_prop(1,2).muav = 0.05;
%  tiss_prop(1,2).gv = 0;
optpos = [0 0 0.1; 1.0 0 0.1; 1.5 0 0.1; 2 0 0.1 ;2.5 0 0.1; 3.0 0 0.1 ];
%optpos = [0 0 0.1;  0.5 0.5 0.1 ];
%optpos = [50 50 2;  75 50 2];
tic
mcxyz_ejecutor('sample', 60,[0 5e-9 5e-9] ,[100 100 100], optpos, tiss_prop, 1);
tiempo_ejec=toc;
lambda=852e-7;
% tauini=1e-7;
% taufinal=1e-2;
% tau=tauini:1e-6:taufinal;
% Db=[1e-8 1.5e-8];
%tic
%[nPhotons,g1,mc]=correlacion_p_mn('sample',lambda,Db,tau,1,1);
mc= mcxyz_read_inp(['sample' num2str(1)]);
s=zeros(3,mc.num_det(1));
for k=1:mc.num_det(1);
    s(:,k)=mc.srcpos(1,:);
end
d=s-mc.detpos';
df=sqrt(sum(d.*d)); %Distance vector source-detector 

%Teorical Auto-correlation function for a semi-infinitum medium
 %g1tcsminfinito=cor_teorica_semiinfinito(mc.tiss_prop(1),mc.tiss_prop(2),mc.tiss_prop(3),1,df,0,lambda,Db,tau);
  
%  for i=1:mc.num_det(1);
% h= figure;
%   semilogx(tau,g1(i,:),'bo',tau,g1tcsminfinito(i,:),'-r')
%   title(['Distancia Fonte-detector=',num2str(df(i)),'cm'],'fontsize',15);
%   xlabel('tau(s)','fontsize',15); ylabel('g1','fontsize',15);
%   axis ([tauini taufinal  0 1 ])
%  h1= legend('montecalo','Semiinfinito');
%  set(h1,'FontSize',14);
%   saveas(h,['correlacionsemiimnfi_' num2str(df(i)) '.png'])
%    saveas(h,['correlacionsemiimnfi_' num2str(df(i)) '.fig'])
%  end
% %time2=toc;
% 
save('sample.mat');