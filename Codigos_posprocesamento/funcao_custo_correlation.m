function chi_correlation= funcao_custo_correlation(x,dados,tau,data_his,ua,Ndet,n,lambda,tiss_num,m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%m Detector
g1=correlacion_mc(x,tau,data_his,0,ua,Ndet,n,lambda,tiss_num);
chi_correlation=zeros(length(g1(m,:)),1);
chi_correlation(1:length(g1(m,:)))=(g1(m,:)-dados)';

end

