
function [momento_T_org pathlength_T_org]=organizador_vetores(momento,pathlenght)

% load('sample.mat');
% %Db=[1e-5 1e-8];
% [pathlenght,momento, nPhoton] = readMCHis_L_M('sample',[mc.tiss_prop(1,3) mc.tiss_prop(2,3)],1 );
% momento=squeeze(momento);
% pathlenght=squeeze(pathlenght);
mom_p=momento(:,1);
%mom_p_copia=momento(:,1);
path_p=pathlenght(:,1);
moment_p_org=sort(mom_p);

i=0;
for j=1:length(mom_p)
  % x=find(moment_p_org==mom_p(j));
    x=find(mom_p==moment_p_org(j));
   k=length(x);
   w=1;
   l=1;
%    if (k==2)
%        break;
%    else
   while w<=k
       i=i+1;
       lst(i)=x(l);
       l=l+1;
       w=w+1;
%        if (l>2)
%            mom_p(x(l-1))=0;
%        end
%            
 %  end
   end
end

for g=1:length(lst)
xx=   find(lst==lst(g));
kk=length(xx);
if kk>1
   e=2;
   while e<=kk
       
    lst(xx(e))=0;
    e=e+1;
   end
    
end
end
lstn=lst(find(lst));
size_m=size(momento);

momento_T_org=zeros(size_m(1),size_m(2));
pathlength_T_org=zeros(size_m(1),size_m(2));
%momento_T_org(:,1)=moment_p_org;

for k=1:size_m(2)
   temp_M= momento(:,k);
   temp_L=pathlenght(:,k);
   
  momento_T_org(:,k)= temp_M(lstn); 
  pathlength_T_org(:,k)=temp_L(lstn);
end
return;
