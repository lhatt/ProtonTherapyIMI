function lhd=likelihoodI(DPa,DP)
%Calculate the likelihood function for a particle

Cd=1e-4*eye(length(DP));% noise matrix corresponding to measurement error
lhd=exp(-.5*(DP-DPa)'*inv(Cd)*(DP-DPa));% likelihood function
