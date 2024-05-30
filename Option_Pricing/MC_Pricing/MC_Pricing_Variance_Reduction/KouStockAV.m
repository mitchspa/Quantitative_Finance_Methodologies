function [S, SAV] = KouStockAV(Nsim, T, params, M, S0, r)
% S = S0 e^{drift+sigmaWt + Jumps}
dt = T/M;
sigma = params(1);
lambda = params(2);
p = params(3);
lamp = params(4);
lamm = params(5);

X = zeros(Nsim, M+1);
XAV = zeros(Nsim, M+1);
psi = @(u) -sigma^2/2*u^2+1i*u*lambda*(p/(lamp-1i*u)-(1-p)/(lamm+1i*u));
drift = r-psi(-1i);

Ndt = icdf('Poisson',rand(Nsim, M),lambda*dt);

Z = randn(Nsim,M);

for i=1:M
    X(:,i+1)=X(:,i)+drift*dt+sqrt(dt)*sigma*Z(:,i);
    XAV(:,i+1)=XAV(:,i)+drift*dt-sqrt(dt)*sigma*Z(:,i);
    for j=1:Nsim

        for jj=1:Ndt(j,i)

            if rand<p %positive
                Y = icdf('Exponential',rand,1/lamp);
                YAV = icdf('Exponential',1-rand,1/lamp);
            else
                Y = -icdf('Exponential',rand,1/lamm);
                YAV = -icdf('Exponential',1-rand,1/lamm);
            
            end

            X(j,i+1)=X(j,i+1)+Y;
            XAV(j,i+1)=XAV(j,i+1)+YAV;
        end
        
    end
end

S = S0*exp(X);
SAV = S0*exp(XAV);