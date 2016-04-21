


%clear all
close all
%clc
format long

% filename='optimizacion_x0=30mm';
% mkdir (filename);
% copyfile('*.m',filename);
% copyfile('mcxyz.exe',filename);
% cd (filename)
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
%rho=1:0.1:3;
%rho=1.0:0.1:1.7;
%                 


load('dados_Dos.mat'),
dados=ap2(1,1:4)/ap2(1,1);


%save('datosexperimentales.mat','dados','rho');
% %%%%%%%%LMsolve%%%%%%%%%%%%%%%%%%%
% x0=30;
% %x0=[0.04 9.5 ]; %4 parametros
% x   = x0;
% x   = x(:);
% lx  = length(x);
% resid=@(x) chi_cuadrado_montecarlo(x,dados,rho); 
% %D   = eval(inp('vektor vah      [D]','1'));
% %ipr = inp('Print control   ipr',0);
% ipr = 1;
% 
% %resid = eval(['@' fun]);    %   function handle
% options = LMFsolve('default');
% % options = LMFsolve...
% %     (options,...
% %     'XTol',1e-6,...
% %     'FTol',1e-6,...
% %     'ScaleD',D,...
% %     'Display',ipr...
% %     );
% 
% 
% options = LMFsolve...
%     (options,...
%     'XTol',1e-1,...
%     'FTol',1e-3,...
%     'Display',ipr,...
%     'Nam', 'opt_3'...
%     );
% 
% %
% [x,S,cnt]=LMFsolve(resid,x,options);

%%%Algoritmo geneticos

% resid=@(x) chi_cuadrado_montecarlo_ga(x,dados,rho,'mc_algoritmo_genetico');
% numberOfVariables = 2;
% [x,fval] = ga(resid,numberOfVariables,[],[],[],[],[0.01 8], [0.1 12]);
% save('mc_genetico.mat')

gaDat3.Objfun='chi_cuadrado_montecarlo_ga';
lb=[0.03 8 0.09 3 ];
ub=[0.08 12 0.35 7];
gaDat3.FieldD=[lb; ub];
% objective function parameters (in this example only used for wasting
% time), see objfun_schwefel_p.m for details.
p.dados=dados;
p.rho=rho;
gaDat3.MAXGEN=20;
gaDat3.ObjfunPar=p;
% Execute GA
gaDat3=ga(gaDat3);
% Result are in
gaDat3.xmin
gaDat3.fxmin
save('genetico.mat')





