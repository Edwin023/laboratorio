%function [x,xs,prob,probs,probt,probts]=pesoscorrelacion(filename,momento,pathlenght,nPhoton,nint,Db,tau,h)
function [x,xs,prob,probs,probt,probts]=pesoscorrelacion(filename,momento,pathlenght,nPhoton,nint,Db,h)
 %[pathlenght,momento, nPhoton] = readMCHis_L_M(filename, mua,h);
mc= mcxyz_read_inp([filename num2str(h)]);
X=zeros(nPhoton,1 );

for i=1:nPhoton
    tem=0;
    for j=1:mc.tiss_num
        tem=tem+6*Db(j)*momento(i,j); 
    
    end
    X(i)=tem; 
end
% 
%  clear
% [pathlenght,momento, nPhoton] = readMCHis_L_M('sample', 0.02,1);
 %nint=round(2*sqrt(nPhoton));
 X=sort(X);
%Divide X em intevalos para fazer o histograma 
X_mom=linspace(min(X),max(X),nint+1);
%mom=linspace(0,max(momento),250);
prob=zeros(1,(length(X_mom)-1));
x=zeros(1,(length(X_mom)-1));
for i=1:(length(X_mom)-1);
    temp=find(X>=X_mom(i) & X<=X_mom(i+1));
    x(i)=((X_mom(i)+X_mom(i+1))/2);
    %Conta  quantos fotones está no intervalo
    prob(i)=length(temp);
    
end
%Encuentra valores no nulos 
tmp=find(prob);
probn=prob(tmp);
 l=1;
  k=probn(1);
  %Cria vetor que repite cada entrada a qantidade de dita entrada
 probt=zeros(1,k);
for j=2:(length(probn))
   
for i=l:k  
       probt(i)=probn(j-1); 
       
end
     
     l=l+probn(j-1);
     k=k+probn(j);

end
for i=l:k
probt(i)=probn(end);
end

probt=probt/nPhoton;
probt=probt/sum(probt);
% figure
% plot(x,prob,'*')


%Agora o mesmo para os comprimentos de caminho

Xs=zeros(nPhoton,1 );
for i=1:nPhoton
    tem=0;
    for j=1:mc.tiss_num
    tem=tem+6*Db(j)*pathlenght(i,j); 
    end
    Xs(i)=tem;
end



s=linspace(min(Xs),max(Xs),nint+1);
%s=linspace(0,max(pathlenght),250);
xs=zeros(1,length(s)-1);
probs=zeros(1,length(s)-1);
for i=1:(length(s)-1);
    temp=find(Xs>=s(i) & Xs<=s(i+1));
    xs(i)=(s(i)+s(i+1))/2;
    probs(i)=length(temp);
    
end

tmp=find(probs);
probn=probs(tmp);
 l=1;
    k=probn(1);
    probts=zeros(1,k);
for j=2:(length(probn))
   
for i=l:k  
       probts(i)=probn(j-1); 
       
end
     l=l+probn(j-1);
     k=k+probn(j);
     
end
for i=l:k
probts(i)=probn(end);
end

probts=probts/nPhoton;
probts=probts/sum(probts);
% figure
%  plot(x,prob,'*')