clc, clear, close all
addpath(genpath('C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks'));
n = 1000;
gamma = 3;
gamma_c = 1;
d = 25;
d_min = 15;
mu = 0.5;

tic
[A,AA,c,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc

%% COMMUNITY DETECTION
Q = community_louvain(A);
NMI = nmi(c,Q);
fprintf('Numero di comunità rilevate: %d\n', max(Q))
fprintf('Normalized Mutual Information: %4f\n',NMI)
W = trustiness(A);

% Distribuzione dei gradi
d_simm = sum(A,2);

figure(1)
subplot(2,1,1)
histogram(d_simm,'BinLimits',[0 50])
subplot(2,1,2)
histogram(dd,'BinLimits',[0 50])


% Stimare il parametro mu

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
figure(2)
histogram(fractions,20,'BinLimits',[0 1])

%%
clc, clear,close all
n = 1000;
mu = 0.4:0.05:0.9;
mm = length(mu);
gamma = [2 3];
gg = length(gamma);
beta = [1 2];
bb = length(beta);
d_mean = [15, 20, 25];
dd = length(d_mean);
d_min = [9,11,16];

A = zeros(n,n,mm,gg,bb,dd);
c = zeros(n,mm,gg,bb,dd);

for m = 1:mm
    for g = 1:gg
        for b = 1:bb
            for d = 1:dd
                [A(:,:,m,g,b,d),AA,c(:,m,g,b,d),dd] = network_LFR(n,d_mean(d),mu(m),gamma(g), beta(b), d_min(d));
            end
        end
    end
end
%%
nmis = zeros(mm,gg,bb,dd);
for m = 1:mm
    for g = 1:gg
        for b = 1:bb
            for d = 1:dd
                Q = community_louvain(A(:,:,m,g,b,d));
                nmis(m,g,b,d) = nmi(c(:,m,g,b,d),Q);
            end
        end
    end
end



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
p = plot(graph(AA|AA'),'EdgeAlpha',0.1);
for d = 1:n
    highlight(p,d,'MarkerSize',log(Q(d)+1),'NodeColor',color(Q(d),:),'EdgeColor','k','LineWidth',0.1)
    p.NodeLabel = [];
end
% figure(2)
% p1 = plot(digraph(A));
% for i = 1:n
%     highlight(p1,i,'MarkerSize',log(c(i)+1),'NodeColor',color(c(i),:))
%     p1.NodeLabel = [];
% end

