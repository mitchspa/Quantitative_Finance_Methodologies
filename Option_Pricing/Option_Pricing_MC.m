function [Price, CI] = Option_Pricing_MC(Nsim, T, M, par, S0, r, model, K, option_type,putORcall)
S = Underlying_Simulation(Nsim, T, M, par, S0, r, model);
Price = NaN;
CI = NaN;

if strcmp(option_type,'EU')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(S(:,end)-K,0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(K-S(:,end),0);
        [Price,~,CI] = normfit(Payoff);
    end

elseif strcmp(option_type,'Lookback')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(max(S,[],2)-K,0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(K-min(S,[],2),0);
        [Price,~,CI] = normfit(Payoff);
    end

elseif strcmp(option_type,'Asian_A_Fs')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(S(:,end)-mean(S(:,:),2),0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(mean(S(:,:),2)-S(:,end),0);
        [Price,~,CI] = normfit(Payoff);
    end

elseif strcmp(option_type,'Asian_A_Fp')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(mean(S(:,:),2)-K,0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(K - mean(S(:,:),2),0);
        [Price,~,CI] = normfit(Payoff);
    end

elseif strcmp(option_type,'Asian_G_Fs')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(S(:,end)-prod(S(:,:),2).^(1/M),0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(prod(S(:,:),2).^(1/M)-S(:,end),0);
        [Price,~,CI] = normfit(Payoff);
    end

elseif strcmp(option_type,'Asian_G_Fp')
    if strcmp(putORcall,'c')
        Payoff = exp(-r*T)*max(prod(S(:,:),2).^(1/M)-K,0);
        [Price,~,CI] = normfit(Payoff);
    elseif strcmp(putORcall,'p')
        Payoff = exp(-r*T)*max(K - prod(S(:,:),2).^(1/M),0);
        [Price,~,CI] = normfit(Payoff);
    end
else
        disp('Error: Unrecognized option type.');
end
end