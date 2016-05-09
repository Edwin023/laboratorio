clear all
close all
clc
A=load('correlation3.mat');
%Cinco parametros
%x(1)=mua1
%x(2)=mus1
%x(3)=mua2
%x(4)=mus2
%x(5)=l é a espesura primeira camada

%Dados experimentais
% amp=[956.7509 167.1661 50.2202 13.4925];
% phase=[-2.4964 -2.3963 -2.3031 -2.2125];

 %Distâncias
%rho=[1.5 2.01 2.52 2.92];
                %tempos de correlation
                tau=1e-7:1e-6:1e-0;
                %Parametros duas camadas teorica 
                n=1; 
                v = 3*10^10/n; % em cm/s
                w =0; % em 1/s
                %aDb1=0; 
                lambda=750e-7;
                %tau=0;
                rho=1.5;
                %aDb2=0;
                cutoff=500;
                mua1=0.05; % em (cm^-1)
                mus1=10;    % em (cm^-1)
                mua2=0.01; % em (cm^-1)
                mus2=5; % em (cm^-1)
                ell=0.7;   % em (cm)
                
                %x_dados=[mua1 mus1 mua2 mus2 ell ];
% %fluencia                
%  dados=phase_amplitud_duascamanas_multi_dist(x_dados,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);  
%correlation
%dados= phase_amplitud_duascamadasref_multi_dist(x_dados,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
dados=A.g1_0(3,:);
% amp=amp/max(amp);
% phi_r=zeros(1,length(rho)-1);
% for j=1:(length(rho)-1)
%     phi_r(j)=dados(2,(j+1))-dados(2,(j));
% end

% dados=zeros(2,length(rho));
% dados(1,:)=amp;
% dados(2,1:length(rho)-1)=phi_r;

x0=[1.1e-10 0.9e-8];


%%%%%%%%%%%%lsqnonlin%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %fluencia
%  phi_duas_camadas=@(x)phase_amplitud_duascamanas_multi_dist(x,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1)-dados;

%Refractancia
 phi_duas_camadas=@(x)correlationforwardsolver_multi_tau(x,n,mua1,mus1,0,tau,lambda,rho,0,w,ell,mua2,mus2,0,cutoff)-dados;

% %trust-region-reflective
%  opts = optimset( 'Display','off');
%  [x_nonlin_t,resnorm_nonlin_t,residuals_nonlin_t,exitflag_nonlin_t,output_nonlin_t] =...
%      lsqnonlin(phi_duas_camadas,x0,[0.001 1 0.001 0.5 0],[0.25 15 0.25 8 2],opts);
%levenberg-marquardt
opts = optimset( 'Algorithm','levenberg-marquardt','Display','off');
[x_nonlin_m,resnorm_nonlin_m,residuals_nonlin_m,exitflag_nonlin_m,output_nonlin_m] = lsqnonlin(phi_duas_camadas,x0,[],[],opts);

%%%%%%%%%%%%lsqcurvefit%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fluencia
% phi_duas_camadas_curvefit=@(x,rho)phase_amplitud_duascamanas_multi_dist(x,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);

%Reflactacia
phi_duas_camadas_curvefit=@(x,tau)correlationforwardsolver_multi_tau(x,n,mua1,mus1,0,tau,lambda,rho,0,w,ell,mua2,mus2,0,cutoff);

% %trust-region-reflective
%   opts = optimset( 'Display','off');
% [x_curvefit_t,resnorm_curvefit_t,residual_curvefit_t] = ...
%     lsqcurvefit(phi_duas_camadas_curvefit,x0,rho,dados,[0.001 1 0.001 0.5 0 ],[0.25 15 0.25 8 2],opts);
%levenberg-marquardt
opts = optimset( 'Algorithm','levenberg-marquardt','Display','off');
[x_curvefit_m,resnorm_curvefit_m,residual_curvefit_m] = lsqcurvefit(phi_duas_camadas_curvefit,x0,tau,dados,[],[],opts);

% %Fluencia
%  phi_opt_nonlin_t=phase_amplitud_duascamanas_multi_dist(x_nonlin_t,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_nonlin_m=phase_amplitud_duascamanas_multi_dist(x_nonlin_m,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_curvefit_t=phase_amplitud_duascamanas_multi_dist(x_curvefit_t,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_curvefit_m=phase_amplitud_duascamanas_multi_dist(x_curvefit_m,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);

%  %Refractancia
%  phi_opt_nonlin_t=phase_amplitud_duascamadasref_multi_dist(x_nonlin_t,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_nonlin_m=phase_amplitud_duascamadasref_multi_dist(x_nonlin_m,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_curvefit_t=phase_amplitud_duascamadasref_multi_dist(x_curvefit_t,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
%  phi_opt_curvefit_m=phase_amplitud_duascamadasref_multi_dist(x_curvefit_m,n,aDb1,tau,lambda,rho,w,aDb2,cutoff,1);
 