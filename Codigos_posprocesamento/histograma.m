
close all
clear
clc
load('sample.mat');
mc = mcxyz_read_inp(['sample' num2str(1)]);

%Np=mc.NPh;
Nd=mc.num_det;

fid = fopen(['sample' num2str(1) '.his'], 'rb');
        data_his = fread(fid, inf, 'float32');
        nhits = size(data_his, 1) / (2*mc.tiss_num+1);
        data_his = reshape(data_his, [2*mc.tiss_num+1 nhits])';
        fclose(fid);
        Nint=[500 350 200 120 70 20];
for k=1:1;
     
 path1=data_his(find(data_his(:,1)==(k-1)), 2:1+mc.tiss_num);
% path1=path1(find(path1(:,1)),:);
  mom1=data_his(find(data_his(:,1)==(k-1)), 2+mc.tiss_num:end);
  nPhoton=size(path1,1);
%Histograma Edwin
[x,xs,prob,probs,probt,probts]=pesoscorrelacion('sample',mom1,path1, nPhoton,Nint(k),1e-8,1);
figure;
% subplot(1,2,1)
% semilogx(x,prob/nPhoton,'*-');
% xlabel('X'); ylabel('P(X)');
% subplot(1,2,2)
semilogx(xs,probs/nPhoton,'*-');
xlabel('s'); ylabel('P(s)');
%   
% %%Histogema Matlab
% 
% X=zeros(nPhoton,1 );
% 
% for i=1:nPhoton
%     tem=0;
%     for j=1:mc.tiss_num
%         tem=tem+6*Db(j)*mom1(i,j); 
%     
%     end
%     X(i)=tem; 
% end
% [n,prob2]=pesos_momentos( X,Nint(k) );
% 
% % %bins=linspace(min(X),max(X),200);
% % [n,x]=hist(X,3000);
% % dx=x(3)-x(2);
% % nf=n/sum(n*dx);
% % %n=n/trapz(bins,n);
% % %sum(n*dx)
% % figure;
% % %semilogx(x,(n))
% % plot(x,nf)
% % 
% % prob=zeros(1,nPhoton);
% % k=1;
% % l=nf(1);
% % for i=1:(length(nf)-1)
% %    prob(k:l)=nf(i)*dx;
% %    k=n(i)+1;
% %    l=n(i+1);
% %    
% % end
% % 
% %  % x_M(k,:)=x;
% % % xs_M(k,:)=xs;
% % % prob_M(k,:)=prob;
% % % probs_M(k,:)=probs;
% % % probt_M(k,:)=probt;
% % % probts_M(k,:)=probts;
 end 
 
 
 
 
 