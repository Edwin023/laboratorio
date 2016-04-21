
function phi=phi_monte_carlo(filename)
% load('mc.mat')
% ua=x(1,:);
% us=x(2,:);

%filename='mc_one_layer_steady';
% copyfile('*.m',filename);
% copyfile('mcxyz.exe',filename);
% cd (filename)

% time_min=3;
% %times_gates=[0 5.0e-9 5.0e-9];
% dim=[100 100 100];
% for i=1:length(us)
% tiss_prop(1,i).musv = us;
% tiss_prop(1,i).muav = ua;
% end
% 
% tiss_prop(1).gv = 0;
% g=tiss_prop(1).gv;

%optpos = [50 50 10;  55 50 10;  60 50 10;  65 50 10;  70 50 10];
%filename='sample_30';

% optpos(i,:)=[0 0 0.1];
% 
% for i=2:length(rho)+1
%   
%    optpos(i,:)=[rho(i-1) 0 0.1]; 
% end

% src=optpos(1,:);
% det=optpos(2:end,:);
% ndet=size(det,1);
% %mcxyz_ejecutor(filename, time_min, dim, optpos, tiss_prop, 1);
% 
% nx =dim(1);
% ny = dim(2);
% nz = dim(3);

% phi_2pt=fluencia_detector_2pt(filename,ndet,nx,ny,nz,det,[0.1 0.1],1);
% 
% phi=zeros(1,length(rho));
% %phi_2pt_max = max(phi_2pt);
%   for i=1:(length(rho)-1)
%     phi(i)=(phi_2pt(i))/(phi_2pt(i+1));
%   end
%  phi(length(rho))=(phi_2pt(1))/(phi_2pt(length(rho)));

[phi_his, I]=mcxyz_read_his(filename,1);

% phi=zeros(1,length(rho));
% %phi_2pt_max = max(phi_2pt);
%   for i=1:(length(rho)-1)
%     phi(i)=log((phi_his(i)))/log((phi_his(i+1)));
%   end
%  phi(length(rho))=log((phi_his(1)))/log((phi_his(length(rho))));

phi=phi_his/max(phi_his);

