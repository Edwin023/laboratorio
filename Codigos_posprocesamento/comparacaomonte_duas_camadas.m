% Comparacion edwin -dois camadas
%clc
close all
%clear
load('sample.mat');
load('tau.mat')
 
%Db=[1e-8 1e-5];
Db=[1e-8 1e-8];
%tau=1e-7:1e-6:1e-2;
tau=taufit;
tau=tau(1:147);
lambda=852e-7;
%clear g1_0;
[~,g1_0,~]=correlacion_p_mn('sample',lambda,Db,tau,1,1,0);
[~,g1_1,~]=correlacion_p_mn('sample',lambda,Db,tau,1,1,1);
[~,g1_2,~]=correlacion_p_mn('sample',lambda,Db,tau,1,1,2);
% tamano=get(0,'ScreenSize');
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% %g1tcsminfinito1=cor_teorica_semiinfinito(mc.tiss_prop(3),1mc.tiss_prop(1,4).0,mc.tiss_prop(2),1,df,0,lambda,Db,tau);
% 
g1_duas=zeros(6,length(tau));
for k=1:6
    for i=1:length(tau)
    %g1_dois(i) = diffusionforwardsolver(mc.tiss_prop(1,4),0,mc.tiss_prop(1,3)*10,mc.tiss_prop(1,1)*10,mc.tiss_prop(1,2),Db(1)/100,tau(i),lambda/10,df/10,0,0.8,mc.tiss_prop(2,3)*10,mc.tiss_prop(2,1)*10,mc.tiss_prop(2,2),Db(2)/100,500);
    %g1_dois(i) = diffusionforwardsolver(mc.tiss_prop(1,4),0,mc.tiss_prop(1,3)*10,mc.tiss_prop(1,1)*10,mc.tiss_prop(1,2),Db(1)/100,tau(i),lambda/10,df/10,0,0.8,mc.tiss_prop(2,3)*10,mc.tiss_prop(2,1)*10,mc.tiss_prop(2,2),Db(2)/100,500);
    g1_duas(k,i) = diffusionforwardsolver(1,0,0.05,10,0,Db(1),tau(i),lambda,df(k),0,0.8,0.05,10,0,Db(2),500);
    end
end
maximo=max(g1_duas,[],2);
for i=1:6
g1_duas(i,:)=g1_duas(i,:)/maximo(i);
end

%%erro tres algoritmos 



% %plot(log10(tau),g1n133(1,:),'bo',log10(tau),g1tcsminfinito1,'-r')
% colores=['b' 'k' 'g' 'r' 'm' 'y'];
% for k=1:1:ndet
% %h=semilogx((tau),(g1_0(1,:)+g1_1(1,:))/2,'bx',(tau),g1_dois,'-g
% %h=semilogx((tau),g1_0(k,:),[colores(k) 'o'],(tau),g1_dois(k,:),['-' colores(k)]);
% h=semilogx((tau),g1_0(k,:),[colores(k) 'o']);
% %h=semilogx((tau),g1_0(1,:),'bx',(tau),g1_dois,'-g');
% set(h(2),'LineWidth',2);
% set(h(1),'MarkerSize',10);
% hh=[xlabel( '$\tau (s)$' ,'interpreter','latex') ylabel( '$g_1$' ,'interpreter','latex')];
% set(hh,'fontsize',16)
% set(gca, 'fontsize', 15) 
% hold on
% %legend('montecarlo','Duas camadas')
% end


% %Variacion us
% load('sample.mat');
% us=0.9:0.01:1.3;
% g1_dois=zeros(length(us),length(tau));
% g1mont=zeros(length(us),length(tau));
% for i=1:length(us)
%     for j=1:length(tau)
%     g1_dois(i,j) = diffusionforwardsolver(mc.tiss_prop(1,4),0,us(i)*10,mc.tiss_prop(1,1)*10,mc.tiss_prop(1,2),Db(1)/100,tau(j),lambda/10,df/10,0,0.8,mc.tiss_prop(2,3)*10,mc.tiss_prop(2,1)*10,mc.tiss_prop(2,2),Db(2)/100,500);
%     
%     end
%     g1mont(i,:)=g1(:,:);
% end
% cout=1e-2;
% dif=(g1_dois-g1mont).*(g1_dois-g1mont)./(g1_dois);
% %dif2=abs((g1tcsminfinito1-g1mont))./(g1mont);
% tt=find(tau<=cout);
% chi2_us=sum(dif(:,1:tt(end)),2);
% tamano=get(0,'ScreenSize');
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% h=plot(us,chi2_us);
% hh=[ylabel('$\chi^2$' ,'interpreter','latex') xlabel(' $\mu_{s}$','interpreter','latex')];
% set(hh,'fontsize',15)
% %saveas(h,['chi_vs_us_prl=' num2str(l) '_cout=' num2str(cout) '.png'])
