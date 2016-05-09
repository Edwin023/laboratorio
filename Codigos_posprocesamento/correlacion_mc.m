function g1=correlacion_mc(x,tau,data_his,l,ua,Nd,n,lambda,tiss_num)
 
k0=2*pi.*n/lambda;
g1=zeros(Nd,length(tau));
for k=1:Nd;
     
 path1=data_his(find(data_his(:,1)==(k-1)), 2:1+mc.tiss_num);
% path1=path1(find(path1(:,1)),:);
  mom1=data_his(find(data_his(:,1)==(k-1)), 2+mc.tiss_num:end);
 %mom1=mom1(find(mom1(:,1)),:);
 %[mom1 path1]=organizador_vetores(mom1,path1); %%%%%P(Y)%%%%%
 nPhotons=size(path1,1);
 
 %momento_total=sum(mom1,2);
 
% w_m=pesos_dis_y(mom1);   %%%Distribucion Y
 
 
    %for  j=1:length(tau)
        X=zeros(nPhotons,1);
        absorcao=zeros(nPhotons,1);
        %dis_2=zeros(nPhotons,1);
        %%%%%%%%%%%%%%%P(X)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %[~,~,~,~,probt,~]=pesoscorrelacion(filenm,mom1,path1,nPhotons, 200,Db ,h);

      
        
        %%%%%%%%%%%%%P(Y)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %[x,xs,prob,probs,probt,probts]=pesoscorrelacion_y(path1,mom1,nPhotons,200);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       for i=1:nPhotons;
            tem=0;
            tem2=0;
         %   temp_dis=0;
            for w=1:tiss_num
                %%%%%%%%%%%%%%P(X)%%%%%%%%%%%%
             %tem=tem+Db(w)*tau(j)*mom1(i,w);
             %tem=tem+Db(w)*mom2(i,w);
             tem=tem+6*x(w)*mom1(i,w);
             tem2=tem2-ua(w)*path1(i,w);
             %%%%%%%%%%%%%%%%%P(Y))%%%%%%%%%%%%%%%%%
            % temp_dis=temp_dis+6*Db(w)*w_m(i,w);
             
            end
%             %%%%%%%%%%P(X)%%%%%%%%%%%%%%%%%%%
               X(i)=tem;
               absorcao(i)=tem2;
           %%%%%%%%%%%%%%%P(Y)%%%%%%%%%%%%%%%%
            %dis_2(i)=temp_dis*momento_total(i);
       end
        %%%%%%%%%%%%%%P(X)%%%%%%%%%%%%%%%%%%%%%%
         [~,probt]=pesos_momentos( X,10*log10(nPhotons));
        [X,absorcao]=organizador_vetores(X,absorcao);
      %[dis_2,absorcao]=organizador_vetores(dis_2,absorcao);

        for i=1:nPhotons;
          %%%%%%%%%%%%%P(X)%%%%%%%%%%%%%%%%%%%%%%%%
          if (l==0)
        g1(k,:)=g1(k,:)+exp((-k0.^2*X(i)).*tau/3).*exp(absorcao(i)); 
          elseif l==1
        g1(k,:)=g1(k,:)+probt(i)*exp(-k0^2*X(i).*tau/3);%+absorcao(i));
        %g1(k,:)=g1(k,:)+probts(i)*exp((-2*(k0^2)*(path1(i)/ls)*Db).*tau);%-ua*path1(i));
          else
        g1(k,:)=g1(k,:)+probt(i)*exp(-k0^2*X(i).*tau/3+absorcao(i));
        %g1(k,:)=g1(k,:)+probts(i)*exp((-2*(k0^2)*(path1(i)/ls)*Db).*tau);%-ua*path1(i));
          end
         %%%%%%%%%%%%%%%%%%%P(Y)%%%%%%%%%%%%%%%%%%%%%
%          if l==0
%         g1(k,:)=g1(k,:)+exp((-k0.^2*dis_2(i).*tau/3)).*exp(absorcao(i)); 
%         % g1(k,:)=g1(k,:)+probt(i)*exp(-k0^2*dis_2(i).*tau/3);%+absorcao(i));
%        %g1(k,j)=g1(k,j)+probt(i)*exp(-2*k0^2*X(i)*tau(j));%-ua(j).*path1(i,j));
%          elseif l==1
%        g1(k,:)=g1(k,:)+probt(i)*exp(-k0^2*dis_2(i).*tau/3);%+absorcao(i));
%         %g1(k,:)=g1(k,:)+probts(i)*exp((-2*(k0^2)*(path1(i)/ls)*Db).*tau);%-ua*path1(i));
%          end
        end
end    

max1=max(g1,[],2)';
for i=1:Nd;
g1(i,:)=g1(i,:)/max1(i);
end

return