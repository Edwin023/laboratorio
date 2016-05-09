function g1 = correlationforwardsolver_multi_tau(x,n,mua1,mus1,go1,tau,lambda,rho,z,w,ell,mua2,mus2,go2,cutoff)
g1=zeros(1,length(tau));
for i=1:length(tau)
g1(i) = diffusionforwardsolver_mod(x,n,mua1,mus1,go1,tau(i),lambda,rho,z,w,ell,mua2,mus2,go2,cutoff);
end
