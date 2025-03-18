
%% 
clc, clear, close all
addpath(genpath('C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks'));
n = 1000;
gamma = 3;
gamma_c = 2;
d = 12;
d_min = 7;
mu = 0.9;

tic
[A,AA,c,dd,s] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc
%%
Q = community_louvain(A);
NMI = nmi(c,Q);
fprintf('Numero di comunità rilevate: %d\n', max(Q))
fprintf('Normalized Mutual Information: %4f\n',NMI)
tic
W = trustiness(A);
toc

%% Distribuzione dei nodi
d_real = sum(A,2);
d_simm = sum((A+A'),2);
figure(1)

subplot(2,1,1)
histogram(d_simm,'BinLimits',[0 50])
subplot(2,1,2)
histogram(dd,'BinLimits',[0 50])




%%
gamma = [3,2];
beta = [1,2];
d_mean = [15 20 25];
mu = 0.4:0.05:0.9;
n = 1000;

for g = gamma
    for b = beta
        for d = d_mean
            for m = mu
            end
        end
    end
end

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
p = plot(digraph(A));
for i = 1:n
    highlight(p,i,'MarkerSize',log(Q(i)+1),'NodeColor',color(Q(i),:),'EdgeColor','k','LineWidth',0.1)
    p.NodeLabel = [];
end
% figure(2)
% p1 = plot(digraph(A));
% for i = 1:n
%     highlight(p1,i,'MarkerSize',log(c(i)+1),'NodeColor',color(c(i),:))
%     p1.NodeLabel = [];
% end

%%
for i = 1:N

    deg = round( (mu)*dd(c == i));
    degrees = sort(deg,'descend');
    degrees = unique(degrees);

    L = zeros(1,n); 

    for j =  1:( length(degrees)-1 )
        d = degrees(j);
        L = L + deg(deg == d);

        if sum(L > 0) < 2
            continue;
        end

        while min(L) > degrees(j+1)

                first_pick = randi([1, sum(L)]);
                % trovo il nodo a cui corrisponde la prima scelta
                k = find(cumsum(L) >= first_pick, 1);
                L(k) = L(k)-1;
                
                % faccio in modo di scegliere un nodo diverso da quello
                % precedente
                L_mod = L;
                L_mod(k) = 0; 

                second_pick = randi([1, sum(L_mod)]);
                % trovo il nodo a cui corrisponde la seconda scelta
                h = find(cumsum(L_mod) >= second_pick, 1);

                % rimuovo i collegamenti effettuati dalla lista
                L(h) = L(h) - 1;

                % creo il collegmento tra i due nodi
                A(k,h) = A(k,h) + 1;
        end
    end
        L = L + deg(deg == degrees(end));
    while sum(L) > 1
            first_pick = randi([1, sum(L)]);
            % trovo il nodo a cui corrisponde la prima scelta
            k = find(cumsum(L) >= first_pick, 1);
            L(k) = L(k)-1;
            
            % faccio in modo di scegliere un nodo diverso da quello
            % precedente
            L_mod = L;
            L_mod(k) = 0; 

            second_pick = randi([1, sum(L_mod)]);
            % trovo il nodo a cui corrisponde la seconda scelta
            h = find(cumsum(L_mod) >= second_pick, 1);

            % rimuovo i collegamenti effettuati dalla lista
            L(h) = L(h) - 1;

            % creo il collegmento tra i due nodi
            A(k,h) = A(k,h) + 1;
    end

            disp(sum(L))
end
