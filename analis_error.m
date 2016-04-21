close all 
load('genetico.mat')
% fid=fopen('puntos_minimos_genetico.txt','r');
%  min_int=fscanf(fid,'%g %g',[2 inf] );
%  min_int=min_int';
%  fclose(fid);
 min_int=gaDat3.xmingen;
 xreal=[0.05 10];
 error=zeros(1,21);
 for i=1:21
     error(i)=sqrt(((min_int(i,1)-xreal(1))/min_int(i,1))^2+((min_int(i,2)-xreal(2))/min_int(i,2))^2);
 end
 figure
 plot(1:21,error)
 xlabel('Generaçao','fontsize',20);
 ylabel('Erro','fontsize',20);
 
 
min_fun=gaDat3.fxmingen; 
%  fid=fopen('Valor_funciones_genetico.txt','r');
%   min_fun=fscanf(fid,'%g',[1 inf] );
%   min_fun=min_fun';
%   
  figure
 plot(1:21,min_fun)
 xlabel('Generaçao','fontsize',20);
 ylabel('valor funçao de custo','fontsize',20);