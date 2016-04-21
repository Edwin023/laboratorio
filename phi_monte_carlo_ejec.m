
function phi=phi_monte_carlo_ejec(x,rho,filename)
%ua=x(1,:);
%us=x(2,:);
ua=zeros(1,2);
us=zeros(1,2);
ua(1)=x(1);
ua(2)=x(3);
us(1)=x(2);
us(2)=x(4);
%ua=0.05;

%filename='mc_one_layer' ;
% copyfile('*.m',filename);
% copyfile('mcxyz.exe',filename);
% cd (filename)

Nphotons=1000000;
%times_gates=[0 5.0e-9 5.0e-9];
rho1=1:0.1:3;
dim=[90 90 90];
for i=1:length(us)
tiss_prop(1,i).musv = us(i);
tiss_prop(1,i).muav = ua(i);
end

tiss_prop(1,1).gv = 0;
tiss_prop(1,2).gv = 0;
%g=tiss_prop(1).gv;

%optpos = [50 50 10;  55 50 10;  60 50 10;  65 50 10;  70 50 10];
%filename='sample_30';

optpos(1,:)=[0 0 0.1];

for j=2:length(rho1)+1
  
   optpos(j,:)=[rho1(j-1) 0 0.1]; 
end

src=optpos(1,:);
det=optpos(2:end,:);
ndet=size(det,1);
mcxyz_ejecutor(filename, Nphotons,dim, optpos, tiss_prop);

nx =dim(1);
ny = dim(2);
nz = dim(3);

% phi_2pt=fluencia_detector_2pt(filename,ndet,nx,ny,nz,det,[0.1 0.1],1);
% 
% phi=zeros(1,length(rho));
% %phi_2pt_max = max(phi_2pt);
%   for i=1:(length(rho)-1)
%     phi(i)=(phi_2pt(i))/(phi_2pt(i+1));
%   end
%  phi(length(rho))=(phi_2pt(1))/(phi_2pt(length(rho)));
 
 
 [phi_his, ~]=mcxyz_read_his(filename);

% phi=zeros(1,length(rho));
% %phi_2pt_max = max(phi_2pt);
%   for i=1:(length(rho)-1)
%     phi(i)=log((phi_his(i)))/log((phi_his(i+1)));
%   end
%  phi(length(rho))=log((phi_his(1)))/log((phi_his(length(rho))));
list=zeros(1,length(rho));

for kk=1:length(rho)
list(kk)=find(rho1==rho(kk));
end
phiT=phi_his/phi_his(list(1));
phi=phiT(list);
%phi=phi_his/max(phi_his);
 %save('mc.mat');
end
