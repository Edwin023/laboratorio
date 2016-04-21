

%function phi = phase_amplitud_duascamadasref_multi_dist(x,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,tmp)
function min=chi_cuadrado_montecarlo_ga(x,p)
dados=p.dados;
rho=p.rho;
filename=p.filename;
% datos=load('datosexperimentales.mat');
% rho=datos.rho;
%A=load('optimizacion.mat');
phi =phi_monte_carlo_ejec(x,rho,filename);


%dados=datos.dados;

%f=phi-dados;
%f=zeros(length(phi)+1,1);
f=zeros(length(phi),1);
%R=(sqrt(x(1)^2/(0.1^2)+x(2)^2/(15^2)+x(3)^2/(0.3^2)+x(4)^2/(6^2))-1).*(x(1)>0 && x(2)>0 && x(3)>0 && x(4)>0);
f(1:length(phi))=(phi-dados)';
sum=0;
for i=1:length(f)
   sum=sum+f(i)*f(i); 
end
min=sum;
%f(length(phi)+1)=(R>0).*R*1000;

return




