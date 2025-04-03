%% PARAMETRI E GENERAZIONE RETE

clc; clear; close all;

n = 10000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.85;

tic
[A,~,c,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc

%% INDIVIDUAZIONE COMUNITÀ

Q_LFR = community_louvain(A);

NMI_LFR = nmi(c,Q_LFR);

fprintf('Numero di comunità rilevate: %d\n', max(Q_LFR));
fprintf('Normalized Mutual Information: %4f\n',NMI_LFR)

%% ANALISI COEFFICIENTE MU

M = (c == c');

sameCommCounts1 = sum(A .* M, 2);
degrees1 = sum(A, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
avgFraction1 = mean(fractions1);
histogram(fractions1,50,"BinLimits",[0 1]);

%% ANALISI GRADO NODI

subplot(2,1,1)
histogram(dd,'BinLimits',[0,50])
subplot(2,1,2)
histogram(sum(A,2),'BinLimits',[0,50])





