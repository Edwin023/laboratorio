

%clear all
close all
clc

%Cinco parametros
%x(1)=mua1
%x(2)=mus1
%x(3)=mua2
%x(4)=mus2
%x(5)=l é a espesura primeira camada

% % % Dados experimentais
% % amp=[5.41498956880929 3.39981546288411 2.92671993875968 53.8808124199552];
% % phase=[-0.109364512591529 -0.109017628402696 -0.0889245253891109 -0.307306666383337];

 %Distâncias
rho=[1.4 1.9 2.5 2.9];
%rho=1.0:0.1:1.7;
%                 Parametros duas camadas teorica 
%                 n=1.33; 
%                 v = 3*10^10/n; % em cm/s
%                 w = 2*pi*110*10^6; % em 1/s
%                 aDb1=0; 
%                 lambda=750e-7;
%                 tau=0;
%                 aDb2=0;
%                 cutoff=500;
                mua=0.05; % em (cm^-1)
                mus=10;    % em (cm^-1)
%                 %ell=1.0;   % em (cm)
%                 
%                 %x_dados=[mua1 mus1 mua2 mus2 ell 1 ];
%                 
                 %x_dados=[mua1 mus1 mua2 mus2 1]; %cinco parametros
                 %x_dados=[mua1 mus1 mua2 mus2]; %cuatro parametros
% %fluencia                
%  dados=phase_amplitud_duascamanas_multi_dist(x_dados,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);  
%Refractancia
%dados= phase_amplitud_duascamadasref_multi_dist(x_dados,rho);


% %dados= phase_amplitud_duascamadas_cilin_ref_multi_dist(x_dados,rho);
% dados=zeros(1,2*length(rho));
% dados(1,1:length(rho))=amp;
% dados(1,length(rho)+1:2*length(rho))=phase;

dados=fluence_fwdsol_v_dist(mua, mus, g, rho, 1) ;



%save('datosexperimentales.mat','dados','rho');


%x0=[0.06 12 0.12 12 1.2 1.0];

%x0=[0.09 11 0.2 3 1.4]; %cinco parametros%
x0=[0.05 10 ]; %4 parametros
%x0=[0.03,5.04,0.109,8.6];
%x0=[0.217341677815244 8.73440829853353 0.0352796573865290 3.35692120872445];
%x   = inp('Starting point [xo]',x0);
x   = x0;
x   = x(:);
lx  = length(x);
resid=@(x) chi_cuadrado_montecarlo(x,dados,rho); 
%D   = eval(inp('vektor vah      [D]','1'));
%ipr = inp('Print control   ipr',0);
ipr = 1;

%resid = eval(['@' fun]);    %   function handle
options = LMFsolve('default');
% options = LMFsolve...
%     (options,...
%     'XTol',1e-6,...
%     'FTol',1e-6,...
%     'ScaleD',D,...
%     'Display',ipr...
%     );


options = LMFsolve...
    (options,...
    'XTol',1e-6,...
    'FTol',1e-6,...
    'Display',ipr...
    );

%
[x,S,cnt]=LMFsolve(resid,x,options);


%save('resultados1.mat')