clear;

%% Parameter Initialization

paramsM = [0.7, 10, -0.01, 0.4];
paramsK = [0.7, 10, 0.5, 10, 10];
Nsim = 1e4; S0 = 1;
T = 1; M = round(52*T); r=0.03; K = 1.2;

%% Option Pricing

[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsM,S0,r, 'merton' ,K, 'Lookback','c');
plot(1,Price,'x',MarkerSize=10)
hold on;
[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsM,S0,r, 'merton', K, 'Asian_A_Fs','c');
plot(1,Price,'x',MarkerSize=10)
[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsM,S0,r, 'merton' ,K, 'Lookback','c');
plot(1,Price,'x',MarkerSize=10)

[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsK,S0,r, 'kou' ,K, 'Lookback','c');
plot(-1,Price,'x',MarkerSize=10)
[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsK,S0,r, 'kou', K, 'Asian_A_Fs','c');
plot(-1,Price,'x',MarkerSize=10)
[Price, CI] = Option_Pricing_MC(Nsim,T,M,paramsK,S0,r, 'kou' ,K, 'Lookback','c');
plot(-1,Price,'x',MarkerSize=10)

legend('EU Call, Merton','Lookback Call, Merton','Asian Floating Strike Call, Merton', ...
    'EU Call, Kou','Lookback Call, Kou','Asian Floating Strike Call, Kou', ...
    'FontSize', 10, 'Location', 'northwest');
title('Warning: random model parameters!')
xlim([-2,2])
ylim([0,1.5])
hold off;
