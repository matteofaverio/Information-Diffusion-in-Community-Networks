clc; clear; close all;

% --- PARAMETRI DEL NETWORK ---
n = 5000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.9;

% --- GENERAZIONE DEL NETWORK ---
tic
[A, c, L, dd] = network_LFRR(n, d, mu, gamma, gamma_c, d_min);
toc
% --- INDIVIDUAZIONE DELLE COMUNITÀ ---(
Q = community_louvain(A);
num_communities = max(Q);
fprintf('Numero di comunità rilevate: %d\n', num_communities)

%%

tic
[AA,err_agg] = rewiring_new(A,c,mu,20000,0.01);
toc

Q1 = community_louvain(AA);
num_communities1 = max(Q1);
fprintf('Numero di comunità rilevate: %d\n', num_communities1);
%%
tic
[AAA,err_agg] = rewiring_new(AA,c,mu,20000,0.01);
toc

Q2 = community_louvain(AAA);
num_communities2 = max(Q2);
fprintf('Numero di comunità rilevate: %d\n', num_communities2);

%%
M = (c == c');

sameCommCounts1 = sum(A .* M, 2);
degrees1 = sum(A, 2);
fractions1 = sameCommCounts1 ./ degrees1;
fractions1(degrees1 == 0) = 0;
avgFraction1 = mean(fractions1);
subplot(2,1,1);
histogram(fractions1,40,"BinLimits",[0 1]);

sameCommCounts = sum(AA .* M, 2);
degrees = sum(AA, 2);
fractions = sameCommCounts ./ degrees;
fractions(degrees == 0) = 0;
avgFraction = mean(fractions);
subplot(2,1,2);
histogram(fractions,40,"BinLimits",[0 1]);



%%

subplot(2,1,1)
histogram(sum(A,2),'BinLimits',[0,50])
subplot(2,1,2)
histogram(sum(AA,2),'BinLimits',[0,50])

%%

it = size(err_agg,2);
for i = 1:it
    errors(i) = mean(err_agg(:,i));
end

plot(errors);
grid on;

%%

randi