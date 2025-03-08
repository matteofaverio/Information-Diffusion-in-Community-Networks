clc; clear; close all;

n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.7;

tic
[B,~,c1,homeless,L1,dd1] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc
fprintf('\n')

tic
[A, c, L, dd] = network_LFRR(n, d, mu, gamma, gamma_c, d_min);
toc
fprintf('\n')

tic
[AA,err_agg,var_agg] = rewiring(A,c,mu,10000);
toc
fprintf('\n')

%%

Q_LFR_rewired = community_louvain(AA);
Q_LFR = community_louvain(B);

NMI_LFR_rewired = nmi(c,Q_LFR_rewired);
NMI_LFR = nmi(c1,Q_LFR);

fprintf('Numero di comunità rilevate con LFR: %d\n', max(Q_LFR));
fprintf('Normalized Mutual Information: %4f\n',NMI_LFR)
fprintf('Numero di comunità rilevate con LFRR+rewiring: %d\n', max(Q_LFR_rewired));
fprintf('Normalized Mutual Information: %4f\n',NMI_LFR_rewired)



%%
M = (c == c');
M1 = (c1 == c1');

sameCommCounts1 = sum(AA .* M, 2);
degrees1 = sum(AA, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
avgFraction1 = mean(fractions1);
subplot(2,1,1);
histogram(fractions1,25,"BinLimits",[0 1]);

sameCommCounts2 = sum(B .* M1, 2);
degrees2 = sum(B, 2);
fractions2 = sameCommCounts2 ./ degrees2;
fractions2(degrees2 == 0) = 0;
avgFraction2 = mean(fractions2);
subplot(2,1,2);
histogram(fractions2,25,"BinLimits",[0 1]);

%%

subplot(2,1,1)
histogram(sum(A,2),'BinLimits',[0,50])
subplot(2,1,2)
histogram(sum(AA,2),'BinLimits',[0,50])

%%

plot(err_agg);
hold on;
grid on;
figure;
plot(var_agg);
grid on;

