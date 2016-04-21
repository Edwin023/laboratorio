function chrom=crtrp_prueba(Nind,FieldDR)
% A random real value matrix is created coerced by upper and 
% lower bounds

Nvar = size(FieldDR,2);
aux = rand(Nind,Nvar);
m=[-1 1]*FieldDR;
ublb=ones(Nind,1)*m;
lb=ones(Nind,1)*FieldDR(1,:);
chrom=ublb.*aux+lb;