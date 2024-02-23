function S = Underlying_Simulation(Nsim, T, M, params, S0, r, model)
% S = S0 e^{drift+sigmaWt+Jumps}

X = zeros(Nsim, M+1);
dt = T/M;

flag = 0;
if strcmp(model,'gbm')
    flag = 1;
elseif strcmp(model,'merton')
    flag = 2;
elseif strcmp(model,'kou')
    flag = 3;
else
    disp('Error in model specification. Use gbm, merton or kou.')
end

switch flag
    case 1
        sigma = params(1);

        psi = @(u) -sigma^2/2*u^2;
    case 2
        sigma = params(1);
        lambda = params(2);
        muJ = params(3);
        sigmaJ = params(4);

        psi = @(u) -sigma^2/2*u^2+lambda*(exp(-sigmaJ^2*u^2/2+1i*muJ*u)-1);
    case 3
        sigma = params(1);
        lambda = params(2);
        p = params(3);
        lamp = params(4);
        lamm = params(5);

        psi = @(u) -sigma^2/2*u^2+1i*u*lambda*(p/(lamp-1i*u)-(1-p)/(lamm+1i*u));
end

drift = r-psi(-1i);

switch flag
    case 1
        
        Z = randn(Nsim,M);
        for i=1:M
            X(:,i+1)=X(:,i)+drift*dt+sqrt(dt)*sigma*Z(:,i);
        end
    case 2
        
        Ndt = icdf('Poisson',rand(Nsim, M),lambda*dt);
        Z = randn(Nsim,M);
        for i=1:M
            X(:,i+1)=X(:,i)+drift*dt+sqrt(dt)*sigma*Z(:,i);
        
            for j=1:Nsim
        
                if Ndt(j,i)>0
        
                    Y = sum(muJ+sigmaJ*randn(Ndt(j,i),1));
                    X(j,i+1)=X(j,i+1)+Y;
                end
            end
        end
    case 3

        Ndt = icdf('Poisson',rand(Nsim, M),lambda*dt);
        Z = randn(Nsim,M);
        for i=1:M
            X(:,i+1)=X(:,i)+drift*dt+sqrt(dt)*sigma*Z(:,i);
        
            for j=1:Nsim
        
                for jj=1:Ndt(j,i)
        
                    if rand<p %positive
                        Y = icdf('Exponential',rand,1/lamp);
                    else
                        Y = -icdf('Exponential',rand,1/lamm);
                    
                    end
        
                    X(j,i+1)=X(j,i+1)+Y;
                end
        
            end
        end

end

S = S0*exp(X);