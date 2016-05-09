clc
%clear all
close all

% temp=input('ingrese 0 para semi-infinito distancia ou 1 para semi-infinito temporal:');
% if temp==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEMI-INFINITO Steady
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='prueba';
copyfile('*.m',filename);
 copyfile('mcxyz.exe',filename);
 cd (filename)

time_min=60;
%times_gates=[0 5.0e-9 5.0e-9];
dim=[100 100 100];
tiss_prop(1).musv = 10;
tiss_prop(1).muav = 0.05;
tiss_prop(1).gv = 0;

tiss_prop(1,2).musv = 10;
tiss_prop(1,2).muav = 0.05;
 tiss_prop(1,2).gv = 0;
% us=tiss_prop(1).musv ;
% ua=tiss_prop(1).muav;
% g=tiss_prop(1).gv;

%optpos = [50 50 10;  55 50 10;  60 50 10;  65 50 10;  70 50 10];
filename='sample_30';

optpos=[0 0 0.1; 0.5 0 0.1; 1.0 0 0.1; 1.5 0 0.1; 2.0  0 0.1; 2.5 0 0.1; 3.0  0 0.1];
%optpos=[0 0 0.1; 0.1 0 0.1;0.2 0 0.1;0.3 0 0.1;0.4 0 0.1;0.5 0 0.1;0.6 0 0.1; 0.7 0 0.1; 0.8 0 0.1; 0.9 0 0.1];
% for i=1:1:21;
%   
%    optpos(i,:)=[0.1*(i-1) 0 0.1]; 
% end

src=optpos(1,:);
det=optpos(2:end,:);
ndet=size(det,1);
  r = sqrt((src(1,1)-det(:,1)).^2 + (src(1,2)-det(:,2)).^2 + (src(1,3)-det(:,3)).^2)';
  %  v = c/n;
  tic;
mcxyz_ejecutor(filename, time_min, dim, optpos, tiss_prop, 1);
time_ejec=toc;
%mc=MCConfig([filename num2str(1)]);
nx =dim(1);
ny = dim(2);
nz = dim(3);

phi_2pt=fluencia_detector_2pt(filename,ndet,nx,ny,det,[0.1 0.1],1);
phi_2pt_max = max(phi_2pt);
phi_2pt_norm = phi_2pt'/phi_2pt_max;
% tic;
% [phi,Nphotons]=tMCimg_read_his(filename,1);
% time_read=toc;
%  phi_his_total=sum(phi);
% phi=phi/phi_his_total;
% phi_his_norm=phi/max(phi);



%fluencia_teorica=fluence_fwdsol_v_dist(ua,us,g,r,1);
%max_fluencia_max=max(fluencia_teorica);
%fluencia_teorica=fluencia_teorica/max_fluencia_max;
save([filename '.mat']);
% h=figure;
% semilogy(r,phi_his_norm,'*',r,fluencia_teorica,'-r');
% xlabel('r(mm)','fontsize',16); ylabel('log (fluencia)','fontsize',16);
% title('Fluencia vs distancia','fontsize',16);
% legend('Montecarlo .his','Semi-infinito','fontsize',15)
% save(filename);
% saveas(h,[filename '.png']);

% h=figure;
% semilogy(r,phi_2pt_norm,'*',r,fluencia_teorica,'-r');
% xlabel('r(mm)','fontsize',16); ylabel('log (fluencia)','fontsize',16);
% title('Fluencia vs distancia','fontsize',16);
% legend('Montecarlo .2pt','Semi-infinito','fontsize',15)
% save(filename);
% saveas(h,[filename '.png']);
% 
% 
% 
% clc
% % close all
% % clear




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %DOS CAPAS-STEADY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename ='Dos_capas_fluence_30d';
%  copyfile('*.m',filename);
%  copyfile('tMCimg.exe',filename);
%  cd (filename)
% 
% NTphotons=1e8;
% times_gates=[0 5.0e-9 5.0e-9];
% dim=[100 100 100];
% tiss_prop(1,1).scattering = 1;
% tiss_prop(1,1).absorption = 0.002;
% tiss_prop(1,1).refraction = 1.4;
% tiss_prop(1,1).anisotropy = 0.8;
% 
% tiss_prop(1,2).scattering = 1.5;
% tiss_prop(1,2).absorption = 0.005;
% tiss_prop(1,2).refraction = 1.4;
% tiss_prop(1,2).anisotropy = 0.8;
% 
% %optpos = [50 50 10;    50 60 10;  50 65 10;  50 70 10; 50 75 10; 50 80 10];
%  for i=1:30;
%     optpos(i,:)=[50 49+i 10]; 
%  end
% 
% src=optpos(1,:);
% det=optpos(2:end,:);
% r = sqrt((src(1,1)-det(:,1)).^2 + (src(1,2)-det(:,2)).^2 + (src(1,3)-det(:,3)).^2)';
% 
% tic
% tMCimg_ejecutor(filename, NTphotons, times_gates, dim, optpos, tiss_prop, 1,'windows');
% tiemeje=toc;
% tic
% [phi,Nphotons]=tMCimg_read_his(filename,1);
% timeread=toc;
%  phi_his_total=sum(phi);
% phi=phi/phi_his_total;
% phi_norm=phi/max(phi);
% file=[filename '.mat'];
% save(file);
% 
% h=figure;
%  semilogy(r,phi_norm,'*');
%  xlabel('t(s)'); ylabel('log (fluencia)');
% file=[filename '.png'];
% saveas(h,file);
% 
% clear all
% clc
% close all
% 
% 
% cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %SIMULACION POLAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NTphotons=1e4;
% % times_gates=[0 5.0e-9 5.0e-9];
% % dim=[100 100 100];
% % tiss_prop(1).scattering = 1.010;
% % tiss_prop(1).absorption = 0.0100;
% % tiss_prop(1).refraction = 1;
% % tiss_prop(1).anisotropy = 0.001;
% % 
% % us=tiss_prop(1).scattering ;
% % ua=tiss_prop(1).absorption;
% % n=tiss_prop(1).refraction;
% % g=tiss_prop(1).anisotropy;
% 
% %optpos = [50 50 10;  55 50 10;  60 50 10;  65 50 10;  70 50 10];
% % filename='sample_polar_30des';
% % mkdir (filename);
% % copyfile('*.m',filename);
% % cd (filename)
% 
% 
% for kkk=1:5
% for kk=1:1
%  %   cd (filename)
%   %  filenm=[filename num2str(kkk)];
%     
% filename=['mc_polar_des_g0_7p_pho_d=' num2str(5*kk) 'mm' num2str(kkk) ];
% mkdir (filename);
% copyfile('*.m',filename);
% copyfile('tMCimg.exe',filename);
% cd (filename)
% 
% NTphotons=1e7;
% times_gates=[0 5.0e-9 5.0e-9];
% dim=[100 100 100];
% tiss_prop(1).scattering = 1.010;
% tiss_prop(1).absorption = 0.0100;
% tiss_prop(1).refraction = 1;
% tiss_prop(1).anisotropy = 0;
% 
% us=tiss_prop(1).scattering ;
% ua=tiss_prop(1).absorption;
% n=tiss_prop(1).refraction;
% g=tiss_prop(1).anisotropy;
% 
% i=2;
% theta=0;
% Ndiv=24;
% optpos=zeros(Ndiv+1,3);
% thetavec(Ndiv)=0;
% optpos(1,:)=[50 50 10];
% while theta<2*pi
%     
%    optpos(i,:)=[50+5*kk*cos(theta) 50+5*kk*sin(theta) 10]; 
%    thetavec(i-1)=theta;
%    theta=theta+2*pi/Ndiv;
%    i=i+1;
% end
% 
% src=optpos(1,:);
% det=optpos(2:end,:);
%   r = sqrt((src(1,1)-det(:,1)).^2 + (src(1,2)-det(:,2)).^2 + (src(1,3)-det(:,3)).^2)';
%   %  v = c/n;
%   tic
% tMCimg_ejecutor(filename, NTphotons, times_gates, dim, optpos, tiss_prop, 1,'windows');
% timeejec=toc;
% tic
% [phi,Nphotons]=tMCimg_read_his(filename,1);
% timeread=toc;
%  phi_his_total=sum(phi);
% phi_TN=phi/phi_his_total;
% phi_norm=phi_TN/max(phi_TN);
% % 
% % 
% % 
% % fluencia_teorica=fluence_fwdsol_v_dist(ua,us,g,r,n);
% % max_fluencia_max=max(fluencia_teorica);
% % fluencia_teorica=fluencia_teorica/max_fluencia_max;
% % 
% %h=figure;
% %polar(thetavec,phi_norm,'-r')
% % title(['fluencia vs Angulo D=' num2str(kk*5) 'mm'],'fontsize',15)
% % saveas(h,['polar' num2str(kk*5) '.png'])
% %file=[filename '.png'];
% %saveas(h,file)
% 
% 
% % saveas(h,['polar' num2str(kk*5) '.png'])
% file=[filename num2str(kkk) '.mat'];
% save(file);
% clear file
% clear
% clc
% cd ..
% end
% end
% % close all


  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Dos capas Temporal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% filename1 ='Dos_capas_fluence_temporal';
%  copyfile('*.m',filename1);
%  copyfile('tMCimg.exe',filename1);
%  cd (filename1)
% 
% NTphotons=1e8;
% times_gates=[0 3.5e-9 1.0e-10];
% dim=[100 100 100];
% tiss_prop(1,1).scattering = 1.1;
% tiss_prop(1,1).absorption = 0.009;
% tiss_prop(1,1).refraction = 1.4;
% tiss_prop(1,1).anisotropy = 0.8;
% 
% tiss_prop(1,2).scattering = 1.3;
% tiss_prop(1,2).absorption = 0.03;
% tiss_prop(1,2).refraction = 1.6;
% tiss_prop(1,2).anisotropy = 0.8;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEMPORAL SEMI-INFINITO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% filename1 ='semi_infinito_fluence_temporal';
%  copyfile('*.m',filename1);
%  copyfile('tMCimg.exe',filename1);
%  cd (filename1)
% 
% 
% NTphotons=1e9;
% times_gates=[0 5.0e-9 1.0e-10];
% dim=[100 100 100];
% tiss_prop(1).scattering = 1.010;
% tiss_prop(1).absorption = 0.0100;
% tiss_prop(1).refraction = 1;
% tiss_prop(1).anisotropy = 0.001;
% 
% us=tiss_prop(1).scattering ;
% ua=tiss_prop(1).absorption;
% n=tiss_prop(1).refraction;
% g=tiss_prop(1).anisotropy;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optpos = [50 50 10;  65 50 10];
% src=optpos(1,:);
% det=optpos(2:end,:);
%  ndet=size(det,1);
%   r = sqrt((src(1,1)-det(:,1)).^2 + (src(1,2)-det(:,2)).^2 + (src(1,3)-det(:,3)).^2)';
%   t=times_gates(3):times_gates(3):times_gates(2);
% 
% tic
% tMCimg_ejecutor(filename1, NTphotons, times_gates, dim, optpos, tiss_prop, 1,'windows');
% timeeje=toc;
% mc=MCConfig([filename1 num2str(1)]);
% nx = mc.ximax - mc.ximin + 1;
% ny = mc.yimax - mc.yimin + 1;
% nz = mc.zimax - mc.zimin + 1;
% nt = mc.nTstep;
% phi_2pt=fluencia_detector_2pt(filename1,ndet,nx,ny,nz,nt,det,1);
% 
% phimean_2pt = mean(phi_2pt,1);
% phi_2pt_max = max(phimean_2pt);
% phi_2pt_norm = phimean_2pt'/phi_2pt_max;
% tic;
% [phi,Nphotons]=tMCimg_read_his(filename1,1);
% time_read=toc;
%  phi_his_total=sum(phi);
% phi=phi/phi_his_total;
% phi_his_norm=phi/max(phi);
% 
% 
% phi_temp = fluence_fwdsol_v_time(ua, us, g, r, t, n);
% phi_temp=phi_temp/max(phi_temp);
% 
% h=figure;
%  semilogy(t,phi_his_norm,'*',t,phi_temp,'-r');
%   xlabel('t(s)','fontsize',16); ylabel('log (fluencia)','fontsize',16);
% title('Fluencia vs tempo','fontsize',16)
% legend('Montecarlo','Semiinfinito','fontsize',15)
%  saveas(h,[filename1 '_his''.png']);
%  
%  
% h=figure;
%  semilogy(t,phi_2pt_norm,'*',t,phi_temp,'-r');
%   xlabel('t(s)','fontsize',16); ylabel('log (fluencia)','fontsize',16);
% title('Fluencia vs tempo','fontsize',16)
% legend('Montecarlo .2pt','Semiinfinito','fontsize',15)
%  saveas(h,[filename1 '_2pt''.png']);
% 
% 
% file=[filename1 '.mat'];
% save(file);
% 



% % h=figure;
% %  semilogy(t,phi_norm,'*');
% %  xlabel('t(s)','fontsize',16); ylabel('log (fluencia)','fontsize'.16);
% %title('Fluencia vs tempo','fontsize',16)
% %  saveas(h,file);
% % 
% % else 
% %     disp('erro as opcoes era 0 um 1')
% % end

