%mc=MCConfig(filename);

%%Parametros
n=1;
Reef=0;
mua1=0.05;
mus1=10;
go1=0;
Db1=1e-8;
rho=1;
lambda=852e-7;
load tau.mat;


mc = mcxyz_read_inp([filenm num2str(h)]);
%Np=mc.NPh;
Nd=mc.num_det;
%n=mc.tiss_prop(:,4); %reflatiction index
%clight0 = 2.99792458e11; %light speed (mm/s)
k0=2*pi.*n/lambda;
us=mc.tiss_prop(:,2)';
g=mc.tiss_prop(:,3)';
ls=1./((1-g).*us);
ua=mc.tiss_prop(:,1)';
%tau=[1e-5:1e-4:1];
 fid = fopen([filenm num2str(h)2 '.his'], 'rb');
        data_his = fread(fid, inf, 'float32');
        nhits = size(data_his, 1) / (2*mc.tiss_num+1);
        data_his = reshape(data_his, [2*mc.tiss_num+1 nhits])';
        fclose(fid);
        
        
  g1_semi = cor_teorica_semiinfinito(x,mua1,mus1,go1,1,rho,0,lambda,taufit);
  
x0=1.2e-8; 
x   = x0;
x   = x(:);
lx  = length(x);
resid=@(x) funcao_custo_correlation(x,g1_semi,taufit,data_his,1 ); 
%D   = eval(inp('vektor vah      [D]','1'));
%ipr = inp('Print control   ipr',0);
ipr = 1;

%resid = eval(['@' fun]);    %   function handle
options = LMFsolve_simples('default');
% options = LMFsolve...
%     (options,...
%     'XTol',1e-6,...
%     'FTol',1e-6,...
%     'ScaleD',D,...
%     'Display',ipr...
%     );


options = LMFsolve...
    (options,...
    'XTol',1e-12,...
    'FTol',1e-12,...
    'Display',ipr...
    );

%
[x,S,cnt]=LMFsolve(resid,x,options);
