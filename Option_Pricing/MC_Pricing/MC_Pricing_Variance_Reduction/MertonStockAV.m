function [S, SAV] = MertonStockAV(Nsim, T, params, M, S0, r)
% S = S0 e^{drift+sigmaWt + Jumps}
dt = T/M;
sigma = params(1);
lambda = params(2);
muJ = params(3);
sigmaJ = params(4);

X = zeros(Nsim, M+1);
XAV = zeros(Nsim, M+1);
psi = @(u) -sigma^2/2*u^2+lambda*(exp(-sigmaJ^2*u^2/2+1i*muJ*u)-1);
drift = r-psi(-1i);

Ndt = icdf('Poisson',rand(Nsim, M),lambda*dt);

Z = randn(Nsim,M);

for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sqrt(dt)*sigma*Z(:,i);
    XAV(:,i+1)=XAV(:,i)+drift*dt-sqrt(dt)*sigma*Z(:,i);
    for j=1:Nsim

        if Ndt(j,i)>0

            Y = sum(muJ+sigmaJ*randn(Ndt(j,i),1));
            YAV = sum(muJ-sigmaJ*randn(Ndt(j,i),1));

            X(j,i+1)=X(j,i+1)+Y;
            XAV(j,i+1)=XAV(j,i+1)+YAV;
        end
    end
end

S = S0*exp(X);
SAV = S0*exp(XAV);
