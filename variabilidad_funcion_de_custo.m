
N=10:10:100;
for j=1:length(N)
filename=['mc_variabilidad_funcion_de_custo' num2str(j)];
copyfile('*.m',filename);
copyfile('mcxyz.exe',filename);
cd (filename)

x=[0.05 10];
min=zeros(1,N(j));
p.rho=1:0.1:3;
dados1=fluence_fwdsol_v_dist(x(1), x(2), 0, rho, 1) ;
p.dados=dados1/max(dados1);
for i=1:N(j)
p.filename=['mc_varifuncusto' num2str(i)];
min(i)=chi_cuadrado_montecarlo_ga(x,p);
end
save(['varabilidad_funcion_de_custo' num2str(j) '.mat'])
cd ..
end