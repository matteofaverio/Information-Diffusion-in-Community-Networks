%% PARAMETRI E GENERAZIONE RETE

clear 


gamma = 3;
gamma_c = 3;
d = 10;
d_min = 7; 

mu_targets = 0.01:0.03:0.99;    % ipotetici valori target per mu_reale
n_values   = 1000:100:2000;  % esempio di numero di nodi

results = struct('mu_target', [], 'n', [], 'A_opt', [], 'k_opt', []);

index = 1;

for n = n_values

    for mu_target = mu_targets

        % Definisci la funzione obiettivo che ora dipende anche da n:
        objective = @(params) objectiveFun(n,d,mu_target,gamma, gamma_c, d_min,params(1),params(2));
        initial_guess = [1, 1];  % ipotesi iniziale per [A, k]
        
        % Esegui l'ottimizzazione (senza vincoli, ad esempio)
        lb = [0, 0];   % limiti inferiori
        ub = [];       % se non ci sono limiti superiori, puoi lasciarlo vuoto o impostare [Inf, Inf]

        %options = optimoptions('fmincon', 'Display', 'iter');  % per vedere l'andamento dell'ottimizzazione
        opt_params = fmincon(objective, initial_guess, [], [], [], [], lb, ub, [], []);
        
        % Salva i risultati
        results(index).mu_target = mu_target;
        results(index).n = n;
        results(index).A_opt = opt_params(1);
        results(index).k_opt = opt_params(2);
        index = index + 1;

    end
end

%%

n = 1000;
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7; 

t = 99;
tests = 5;
tr = 100;

mu = linspace(0.01,0.99,t);
% AVG_MU = zeros(t,tests);
parameters = zeros(t,3);
tic

for i = 1:t

    A = 1;
    k = 1;

    mu_test = zeros(tests,1);
    for j = 1:tr

        
        for s = 1:tests
    
            %tic
            [A,~,c,dd] = LFR2(n,d,mu(i),gamma, gamma_c, d_min, A, k); 
            %toc
            
            M = (c == c');
            sameCommCounts1 = sum(A .* M, 2);
            degrees1 = sum(A, 2);
            fractions1 = sameCommCounts1 ./ degrees1;
            fractions1(degrees1 == 0) = 0;
            mu_test(s) = mean(fractions1);
   
        end

    end
    parameters(i,:) = [mu(i) A k];
    toc
    fprintf('mu = %d\n\n', mu(i));
    tic
end
toc

%%

nonzero_count = sum(AVG_MU ~= 0, 2); 
row_sum = sum(AVG_MU, 2);  
media_righe = row_sum ./ nonzero_count;
nonzero_count1 = sum(AVG_MU1 ~= 0, 2); 
row_sum1 = sum(AVG_MU1, 2);  
media_righe1 = row_sum1 ./ nonzero_count1;

x = mu;
y = mu'-media_righe;

coeff = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);  % crea una griglia fine
y_fit = polyval(coeff, x_fit);
plot(x, y)%, x_fit, y_fit, '-')


%%

err = mu'-media_righe;
err_mediato = zeros(t,1);
h = 6;

for i = 1:t
        idx_start = max(1, i - h);
        idx_end   = min(t, i + h);
        err_mediato(i) = mean(err(idx_start:idx_end));
end
plot(mu,err_mediato);


%% ANALISI GRADO NODI

subplot(2,1,1)
histogram(dd,'BinLimits',[0,50])
subplot(2,1,2)
histogram(sum(A,2),'BinLimits',[0,50])

%% INDIVIDUAZIONE COMUNITÀ

Q_LFR = community_louvain(A);

NMI_LFR = nmi(c,Q_LFR);

fprintf('Numero di comunità rilevate: %d\n', max(Q_LFR));
fprintf('Normalized Mutual Information: %4f\n',NMI_LFR)

%%

clc; clear; close all;

n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.85;
%%
tic
[A,~,c,dd] = LFR2(n,d,mu,gamma, gamma_c, d_min);
toc

%%

www = [816
   939
   963
   967];

%% 

t = [ 2 3 ; 4 3 ; 3 2]

t(t(:,1) == 3,:) = [];

t


