clear all;


paramsM = [0.7, 10, -0.01, 0.4];
paramsK = [0.7, 10, 0.5, 10, 10];
Nsim = 1e4; S0 = 100;
T = 1; M = round(52*T); r=0.03;

%% Stock Simulation
[SM,SM_AV] = MertonStockAV(Nsim,T,paramsM,M,S0,r);
[SK,SK_AV] = KouStockAV(Nsim,T,paramsK,M,S0,r);

%% Option Pricing
K = 102;

Payoff_M = exp(-r*T)*max(SM(:,end)-K,0);
[PayoffM,~,PayoffMCI] = normfit(Payoff_M);

Payoff_K = exp(-r*T)*max(SK(:,end)-K,0);
[PayoffK,~,PayoffKCI] = normfit(Payoff_K);

C_M = [PayoffMCI(1),PayoffM,PayoffMCI(2)];
C_K = [PayoffKCI(1),PayoffK,PayoffKCI(2)];




Payoff_MAV = exp(-r*T)*max(SM_AV(:,end)-K,0);
[PayoffM,~,PayoffMCI] = normfit((Payoff_MAV+Payoff_M)/2);

Payoff_KAV = exp(-r*T)*max(SK_AV(:,end)-K,0);
[PayoffK,~,PayoffKCI] = normfit((Payoff_KAV+Payoff_K)/2);

C_M_AV = [PayoffMCI(1),PayoffM,PayoffMCI(2)];
C_K_AV = [PayoffKCI(1),PayoffK,PayoffKCI(2)];

disp('Standard MC:')
C_M
C_K
disp('Variance Reduction:')
C_M_AV
C_K_AV
