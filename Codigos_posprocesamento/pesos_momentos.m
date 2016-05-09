function  [n,prob]  = pesos_momentos( X,Nint )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

     
%%Histogema Matlab
% nPhoton=size(mom1,1);
% X=zeros(nPhoton,1 );
% 
% for i=1:nPhoton
%     tem=0;
%     for j=1:tiss_num
%         tem=tem+6*Db(j)*mom1(i,j); 
%     
%     end
%     X(i)=tem; 
% end

%bins=linspace(min(X),max(X),200);
[n,x]=hist(X,Nint);
dx=x(3)-x(2);
nf=n/sum(n*dx);
%n=n/trapz(bins,n);
%sum(n*dx)
% figure;
% %semilogx(x,(n))
% plot(x,nf)

prob=zeros(1,length(X));
k=1;
l=n(1);
%sum=n(1);
for i=1:(length(nf)-1)
   prob(k:l)=nf(i)*dx;
   k=k+n(i);
   %k=n(i)+1;
   l=l+n(i+1);
   
end

end

