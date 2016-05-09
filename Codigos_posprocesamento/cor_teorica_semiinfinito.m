function g1=cor_teorica_semiinfinito(x,ua,us,g,n,d,z,lambda,tau)
us=us*(1-g);
if n==1.333

Reef=0.431;
elseif n==1.4
Reef=0.493;
elseif n==1
Reef=0;
elseif n==1.431
Reef=0.431;
end

z0=1/(us+ua);
zb=((1+Reef)/(1-Reef))*2*z0/3;
r1=sqrt(d.^2+(z-z0)^2);
r2=sqrt(d.^2+(z+z0+2*zb)^2);
k0=2*pi*n/lambda;
%k=sqrt(3*us*ua+6*us^2*k0^2*Db*tau);
k=sqrt(3*(us+ua)*ua+6*us*(ua+us)*k0^2*x*tau);
g1=zeros(length(d),length(tau));
for j=1:length(d)
g1(j,:)=(3*(ua+us)/(4*pi))*(exp(-k*r1(j))/r1(j)-exp(-k*r2(j))/r2(j));
end
max1=max(g1,[],2)';
for i=1:length(d);
g1(i,:)=g1(i,:)/max1(i);
end

return

