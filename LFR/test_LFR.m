clc,clear,close all
gamma = 3;
d_min = 10;
n = 10000;
d = 19;
A = zeros(n);

%  generazione casuale dei gradi dei singoli nodi in modo che seguano una
%  power law distribution:

% calcolo di d_max con fissato
num = @(x) (x.^(2-gamma)-d_min.^(2-gamma))/(2-gamma);
den = @(x) (x.^(1-gamma)-d_min.^(1-gamma))/(1-gamma);
dmean = @(x) num(x)./den(x);

% funzione a cui applicare il metodo di bisezione
f = @(x) dmean(x) - d;

% estremi dell'intervallo su cui applicare la bisezione
lower = d_min + 1e-9;
upper = 100.0;

% sposto l'upper bound in modo da soddisfare le condizioni iniziali del
% metodo di bisezione
while f(lower) * f(upper) > 0
    upper = upper * 2;
    if upper > 1e9
        error(['Non riesco a trovare un intervallo adeguato per d_max. ', ...
               'Prova ad altri parametri o a ridurre la media desiderata.']);
    end
end

d_max1 = fzero(f,[lower, upper]);

u = rand(n,1);
dd = ( d_min^(1-gamma) + u*(d_max1^(1-gamma) - d_min^(1-gamma))).^(1/(1-gamma));
%%
n = 100000;
d_min = 40;
d_max = 1000;
gamma_c = 3;
[S,N] = powerLaw_communities(n,d_min,d_max,gamma_c);
figure;
histogram(S, 'Normalization','pdf', 'EdgeColor','none');
title('Istogramma dei campioni ~ Power Law troncata');
xlabel('x');
ylabel('Densità di probabilità stimata');

%% 
clc, clear, close all
addpath(genpath('C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks'));
n = 1000;
gamma = 3;
gamma_c = 2;
d = 22;
d_min = 15;
mu = 0.8;

tic
[A,AA,c,h,L,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc
Q = community_louvain(A);
NMI = nmi(c,Q);
fprintf('Numero di comunità rilevate: %d\n', max(Q))
fprintf('Normalized Mutual Information: %4f\n',NMI)
tic
W = trustiness(A);
toc


%% Stimare il parametro mu

% maschera booleana, M(i,j) == 1 se i e j fanno parte della stessa comunità
M = (c == c');
% applico la maschera ad A e sommo sulle righe: trovo i collegamenti
% interni di ogni nodo
sameCommCounts = sum(A .* M, 2);
% trovo i collegamenti totali dei nodi 
degrees = sum(A, 2);
fractions = sameCommCounts ./ degrees;
% tolgo i nan
fractions(degrees == 0) = 0;
avgFraction = mean(fractions);
fprintf('mu medio rilevato: %f\n',avgFraction)
histogram(fractions,10)

mu_unitario = find(fractions == 1);
nodistrani = dd(mu_unitario);

%%
color = [0.00, 0.45, 0.70;  0.85, 0.33, 0.10;  0.93, 0.69, 0.13;  0.49, 0.18, 0.56;
    0.47, 0.67, 0.19;  0.30, 0.75, 0.93;  0.64, 0.08, 0.18;  0.30, 0.30, 0.30;
    0.60, 0.60, 0.60;  1.00, 0.00, 0.00;  1.00, 0.50, 0.00;  0.75, 0.75, 0.00;
    0.00, 1.00, 0.00;  0.00, 0.00, 1.00;  0.66, 0.00, 1.00;  0.33, 0.33, 0.00;
    0.33, 0.67, 0.00;  0.33, 0.00, 0.67;  0.67, 0.33, 0.00;  0.67, 0.67, 0.33;
    0.33, 0.67, 0.67;  0.67, 0.33, 0.67;  0.33, 0.33, 0.67;  0.67, 0.00, 0.33;
    0.00, 0.67, 0.33;  0.13, 0.55, 0.13;  0.70, 0.13, 0.13;  0.10, 0.25, 0.40;
    0.98, 0.50, 0.45;  0.88, 0.44, 0.84;  0.19, 0.58, 0.78;  0.85, 0.57, 0.94;
    0.92, 0.78, 0.62;  0.80, 0.36, 0.36;  0.40, 0.50, 0.80;  0.12, 0.63, 0.42;
    0.99, 0.80, 0.20;  0.90, 0.30, 0.20;  0.50, 0.20, 0.50;  0.20, 0.80, 0.60;
    0.60, 0.20, 0.80;  0.90, 0.60, 0.30;  0.10, 0.70, 0.70;  0.70, 0.10, 0.70;
    0.75, 0.25, 0.50;  0.50, 0.75, 0.25;  0.25, 0.50, 0.75;  0.60, 0.40, 0.40;
    0.40, 0.60, 0.40;  0.40, 0.40, 0.60;  0.30, 0.30, 0.30;  0.80, 0.50, 0.20;
    0.90, 0.90, 0.20;  0.20, 0.90, 0.90;  0.90, 0.20, 0.90;  0.30, 0.80, 0.30;
    0.80, 0.30, 0.30;  0.30, 0.30, 0.80;  0.70, 0.50, 0.50;  0.50, 0.70, 0.50;
    0.50, 0.50, 0.70;  0.25, 0.75, 0.25;  0.75, 0.25, 0.25;  0.25, 0.25, 0.75;
    0.80, 0.80, 0.40;  0.40, 0.80, 0.80;  0.80, 0.40, 0.80;  0.20, 0.60, 0.20;
    0.60, 0.20, 0.20;  0.20, 0.20, 0.60;  0.90, 0.30, 0.50;  0.50, 0.90, 0.30;
    0.30, 0.50, 0.90;  0.70, 0.70, 0.30;  0.30, 0.70, 0.70;  0.70, 0.30, 0.70;
    0.10, 0.40, 0.10;  0.40, 0.10, 0.10;  0.10, 0.10, 0.40;  0.95, 0.55, 0.15;
    0.15, 0.95, 0.55;  0.55, 0.15, 0.95;  0.80, 0.60, 0.40;  0.40, 0.80, 0.60;
    0.60, 0.40, 0.80;  0.50, 0.80, 0.50;  0.80, 0.50, 0.80;  0.50, 0.50, 0.80;
    0.95, 0.70, 0.30;  0.30, 0.95, 0.70;  0.70, 0.30, 0.95;  0.60, 0.20, 0.40;
    0.20, 0.60, 0.40;  0.40, 0.20, 0.60;  0.30, 0.90, 0.90;  0.90, 0.30, 0.90;
    0.90, 0.90, 0.30;  0.80, 0.80, 0.80;  0.20, 0.20, 0.20;  0.50, 0.50, 0.50;
];
figure(1)
p = plot(digraph(AA));
for i = 1:n
    highlight(p,i,'MarkerSize',log(Q(i)+1),'NodeColor',color(Q(i),:))
    p.NodeLabel = [];
end
% figure(2)
% p1 = plot(digraph(A));
% for i = 1:n
%     highlight(p1,i,'MarkerSize',log(c(i)+1),'NodeColor',color(c(i),:))
%     p1.NodeLabel = [];
% end

%%
clc, clear, close all
n = 100;
gamma = 3;
gamma_c = 3;
d = 7;
d_min = 4;
mu = 0.8;
[A,AA,c,h,L,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
Q = community_louvain(A);
color = [0.00, 0.45, 0.70;  0.85, 0.33, 0.10;  0.93, 0.69, 0.13;  0.49, 0.18, 0.56;
    0.47, 0.67, 0.19;  0.30, 0.75, 0.93;  0.64, 0.08, 0.18;  0.30, 0.30, 0.30;
    0.60, 0.60, 0.60;  1.00, 0.00, 0.00;  1.00, 0.50, 0.00;  0.75, 0.75, 0.00;
    0.00, 1.00, 0.00;  0.00, 0.00, 1.00;  0.66, 0.00, 1.00;  0.33, 0.33, 0.00;
    0.33, 0.67, 0.00;  0.33, 0.00, 0.67;  0.67, 0.33, 0.00;  0.67, 0.67, 0.33;
    0.33, 0.67, 0.67;  0.67, 0.33, 0.67;  0.33, 0.33, 0.67;  0.67, 0.00, 0.33;
    0.00, 0.67, 0.33;  0.13, 0.55, 0.13;  0.70, 0.13, 0.13;  0.10, 0.25, 0.40;
    0.98, 0.50, 0.45;  0.88, 0.44, 0.84;  0.19, 0.58, 0.78;  0.85, 0.57, 0.94;
    0.92, 0.78, 0.62;  0.80, 0.36, 0.36;  0.40, 0.50, 0.80;  0.12, 0.63, 0.42;
    0.99, 0.80, 0.20;  0.90, 0.30, 0.20;  0.50, 0.20, 0.50;  0.20, 0.80, 0.60;
    0.60, 0.20, 0.80;  0.90, 0.60, 0.30;  0.10, 0.70, 0.70;  0.70, 0.10, 0.70;
    0.75, 0.25, 0.50;  0.50, 0.75, 0.25;  0.25, 0.50, 0.75;  0.60, 0.40, 0.40;
    0.40, 0.60, 0.40;  0.40, 0.40, 0.60;  0.30, 0.30, 0.30;  0.80, 0.50, 0.20;
    0.90, 0.90, 0.20;  0.20, 0.90, 0.90;  0.90, 0.20, 0.90;  0.30, 0.80, 0.30;
    0.80, 0.30, 0.30;  0.30, 0.30, 0.80;  0.70, 0.50, 0.50;  0.50, 0.70, 0.50;
    0.50, 0.50, 0.70;  0.25, 0.75, 0.25;  0.75, 0.25, 0.25;  0.25, 0.25, 0.75;
    0.80, 0.80, 0.40;  0.40, 0.80, 0.80;  0.80, 0.40, 0.80;  0.20, 0.60, 0.20;
    0.60, 0.20, 0.20;  0.20, 0.20, 0.60;  0.90, 0.30, 0.50;  0.50, 0.90, 0.30;
    0.30, 0.50, 0.90;  0.70, 0.70, 0.30;  0.30, 0.70, 0.70;  0.70, 0.30, 0.70;
    0.10, 0.40, 0.10;  0.40, 0.10, 0.10;  0.10, 0.10, 0.40;  0.95, 0.55, 0.15;
    0.15, 0.95, 0.55;  0.55, 0.15, 0.95;  0.80, 0.60, 0.40;  0.40, 0.80, 0.60;
    0.60, 0.40, 0.80;  0.50, 0.80, 0.50;  0.80, 0.50, 0.80;  0.50, 0.50, 0.80;
    0.95, 0.70, 0.30;  0.30, 0.95, 0.70;  0.70, 0.30, 0.95;  0.60, 0.20, 0.40;
    0.20, 0.60, 0.40;  0.40, 0.20, 0.60;  0.30, 0.90, 0.90;  0.90, 0.30, 0.90;
    0.90, 0.90, 0.30;  0.80, 0.80, 0.80;  0.20, 0.20, 0.20;  0.50, 0.50, 0.50;
];
p = plot(digraph(A));
for i = 1:n
    highlight(p,i,'MarkerSize',(c(i)+1),'NodeColor',color(Q(i),:))
    p.NodeLabel = [];
end