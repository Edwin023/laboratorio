% close all
% 
% %Variacion us
% load('sample.mat');
% [nPhotons2,g1,mc2]=correlacion_p_mn('sample',lambda,Db,tau,1,1);
% %juan=load('juan.mat');
% us=0.7:0.01:1;
% g1tcsminfinito1=zeros(length(us),length(tau));
% g1mont=zeros(length(us),length(tau));
% for i=1:length(us)
% g1tcsminfinito1(i,:)=cor_teorica_semiinfinito(mc.tiss_prop(3),us(i),mc.tiss_prop(2),1,df,0,lambda,Db,tau);
% g1mont(i,:)=g1(:,:);
% end
% 
%  dif=(g1tcsminfinito1-g1mont).*(g1tcsminfinito1-g1mont)./(g1tcsminfinito1);
%  dif2=abs((g1tcsminfinito1-g1mont))./(g1mont);
%  cout=1e-3;
%   tt=find(tau<=cout);
%   chi2_us=sum(dif(:,1:tt(end)),2);
% 
%  h=figure;
%  plot(us,chi2_us);
%  xlabel('us');
%  ylabel('chi2');
% %  saveas(h,['chi_vs_l_prz=' num2str(z) '_cout=' num2str(cout) '.png'])

% %  
% %  
% % h= figure;
% %  maximo=max(dif2(:,1:tt(end))*100,[],2);
% %  plot(l,maximo);
% %  xlabel('parametro l');
% %  ylabel('maximorelativo');
% %  saveas(h,['maximo_vs_l_prz=' num2str(z) 'cout=' num2str(cout) '.png'])
% %  
% 
% %  Variacion z
% %load('sample.mat');
% 
% z=0:0.5:10;
% g1tcsminfinito1=zeros(length(z),length(tau));
% g1mont=zeros(length(z),length(tau));
% for i=1:length(z)
% g1tcsminfinito1(i,:)=cor_teorica_semiinfinito(mc.tiss_prop(3),2.2,mc.tiss_prop(2),1,df,z(i),lambda,Db,tau);
% g1mont(i,:)=g1(:,:);
% end
% 
%  dif=(g1tcsminfinito1-g1mont).*(g1tcsminfinito1-g1mont)./(g1tcsminfinito1);
%  %dif2=abs((g1tcsminfinito1-g1mont))./(g1mont);
% 
%  tt=find(tau<=cout);
%  chi2_z=sum(dif(:,1:tt(end)),2);
%  figure;
%  
%  h=plot(z,chi2_z);
%  xlabel('parametro z');
%  ylabel('chi2');
% % saveas(h,['chi_vs_z_prl=' num2str(l) '_cout=' num2str(cout) '.png'])
% 
%  
%  
%  figure;
%  maximo=max(dif2(:,1:tt(end))*100,[],2);
%  plot(z,maximo);
%  xlabel('parametro z');
%  ylabel('maximorelativo');
%  saveas(h,['maximo_vs_z_prl=' num2str(l) '_cout=' num2str(cout) '.png'])
% 
%  
%  
%  
% % Comparacion con n=1 
%  g1tcsminfinito1teorica=cor_teorica_semiinfinito(mc.tiss_prop(3),mc.tiss_prop(1),mc.tiss_prop(2),1,df,1,lambda,Db,tau,1);
%  
%  for i=1:mc.num_det(1);
% h= figure;
%   semilogx(tau,g1(i,:),'bo',tau,g1tcsminfinito1teorica(i,:),'-r')
%   title(['Distancia Fonte-detector=',num2str(df(i)),'mm'],'fontsize',15);
%   xlabel('tau(s)','fontsize',15); ylabel('g1','fontsize',15);
%   axis ([tauini taufinal  0 1 ])
%  h1= legend('montecalo','Semiinfinito');
%  set(h1,'FontSize',14);
% %  saveas(h,['correlacionsemiimnfi_' num2str(df(i)) '.png'])
%  end
 
 
 
% % Comparacion juan-edwin
% figure;
% load('sample.mat');
%  juan=load('juan.mat');
%  [nPhotons2,g1n133,mc2]=correlacion_p_mn('sample',lambda,Db,tau,1,1.333);
% 
%  semilogx(tau,g1n133(1,:),'bo',tau,juan.gtaus,'-go')
%  legend('edwin','juan')
 
 
 
  
%  % Comparacion con n=1.33 
%  load('sample.mat');
%  g1tcsminfinito1teorica=cor_teorica_semiinfinito(mc.tiss_prop(3),mc.tiss_prop(1),mc.tiss_prop(2),1.333,df,1,lambda,Db,tau,1);
%  [nPhotons1,g1133,mc1]=correlacion_p_mn('sample',lambda,Db,tau,1,1.333);
%  for i=1:mc.num_det(1);
% h= figure;
%   semilogx(tau,g1133(i,:),'bo',tau,g1tcsminfinito1teorica(i,:),'-r')
%   title(['Distancia Fonte-detector=',num2str(df(i)),'mm'],'fontsize',15);
%   xlabel('$\tau(s)$','interpreter','latex','fontsize',20); ylabel('g1','fontsize',15);
%   axis ([tauini taufinal  0 1 ])
%  h1= legend('montecalo','Semiinfinito');
%  set(h1,'FontSize',14);
% %  saveas(h,['correlacionsemiimnfi_' num2str(df(i)) '.png'])
%  end


% Comparacion edwin -semi
close all
clear all
 clc
load('sample.mat');
tamano=get(0,'ScreenSize');
for k=1:6;
figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);

%[nPhotons2,g1n133,mc2]=correlacion_p_mn('sample',lambda,Db,tau,1,1);
%g1tcsminfinito1=cor_teorica_semiinfinito(mc.tiss_prop(3),1.0,mc.tiss_prop(2),1,df,0,lambda,Db,tau);
% plot(log10(tau),g1n133(1,:),'bo',log10(tau),g1tcsminfinito1,'-r')
 
 h=plot(log10(tau),g1_1(k,:),'bo',log10(tau),g1tcsminfinito(k,:),'-r');
% h=semilogy(sqrt(tau/1e-6),(g1n133(1,:)),'bo',sqrt(tau/1e-6),(g1tcsminfinito1),'-r');
set(h(2),'LineWidth',2);
set(h(1),'MarkerSize',10);
%axis([0 20 1e-10 1])
hh=[xlabel( '$\tau(s)$' ,'interpreter','latex') ylabel('$g_1$','interpreter','latex')];
title(['Funçao de autocorrelaçao para distancia F-D=' num2str(df(k)) 'cm'],'fontsize',18)
set(hh,'fontsize',16)
set(gca, 'fontsize', 13) 
legend('Montecarlo','semi-infinito')

end
