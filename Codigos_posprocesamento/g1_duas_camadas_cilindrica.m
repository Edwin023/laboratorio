%function phi=phi_duas_camadas_cilindrica(x,n,z,a,rho,w,l2,cutoff)
function phi=g1_duas_camadas_cilindrica(x,n,z,a,rho,w,aDb1,aDb2,tau,lambda,l2,cutoff)

if n==1
    Reff=0;
elseif n==1.33
    Reff=0.431;
   elseif n==1.4
    Reff=0.439;
end

D1=1/(3*(x(1)+x(2)));
zb=2*D1*(1+Reff)/(1-Reff);
aprima=a+zb;

S=besselzero(0,cutoff,1)/aprima;
   % S =gl(1:cutoff,2);
    phi=0;
for i=1:cutoff
   % phi=phi+green_duas_camadas(x,S(i),z,w,n,l2).*besselj(0,rho.*S(i))./(besselj(1,aprima.*S(i)).^2);
    phi=phi+green_duas_camadas(x,S(i),z,w,n,aDb1,aDb2,tau,lambda,l2).*besselj(0,rho.*S(i))./(besselj(1,aprima.*S(i)).^2);
end
phi=phi/(pi*aprima^2);
return


function green=green_duas_camadas(x,S,z,w,n,aDb1,aDb2,tau,lambda,l2)
c = 2.9979e10;
v=c/n;
k0=2*pi/lambda;
z0=1/((x(1)+x(2)));
if n==1
    Reff=0;
elseif n==1.33
    Reff=0.431;
   elseif n==1.4
    Reff=0.439;
end

D1=1/(3*(x(1)+x(2)));
D2=1/(3*(x(3)+x(4)));

zb1=2*D1*(1+Reff)/(1-Reff);
zb2=2*D2*(1+Reff)/(1-Reff);

alfa1=x(1)/D1+S.^2+2*x(2)*k0^2.*tau*aDb1/D1+1i*w/(D1*v);
alfa1sq=sqrt(alfa1);
alfa2=x(3)/D2+S.^2+2*x(4)*k0^2.*tau*aDb2/D2+1i*w/(D2*v);
alfa2sq=sqrt(alfa2);

betha3=sinh(alfa2sq.*(l2+zb2));
gamma3=cosh(alfa2sq.*(l2+zb2));

green=(exp(-alfa1sq*abs(z-z0))-exp(-alfa1sq*(z+z0+2*zb1)))./(2*D1*alfa1sq)+(sinh(alfa1sq*(z0+zb1)).*sinh(alfa1sq*(z+zb1))./(D1*alfa1sq.*exp(alfa1sq*(x(5)+zb1))))...
    .*((D1*alfa1sq.*betha3-D2*alfa2sq.*gamma3))./(D1*alfa1sq.*cosh(alfa1sq*(x(5)+zb1)).*betha3+D2*alfa2sq.*sinh(alfa1sq*(x(5)+zb1)).*gamma3);


% green=(exp(-alfa1sq*abs(z-z0))-exp(-alfa1sq*(z+z0+2*zb1)))./(2*D1*alfa1sq)+(sinh(alfa1sq*(z0+zb1)).*sinh(alfa1sq*(z+zb1))./(D1*alfa1sq.*exp(alfa1sq*(L+zb1))))...
%     .*((D1*alfa1sq.*betha3-D2*alfa2sq.*gamma3))./(D1*alfa1sq.*cosh(alfa1sq*(L+zb1)).*betha3+D2*alfa2sq.*sinh(alfa1sq*(L+zb1)).*gamma3);

return