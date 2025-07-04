
%% Generazione Rete
n = 2000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;

t = 1;
TEST = cell(t,t);
mu = 0.9; %linspace(0.3,0.999,t);

flag = true;
tic
while flag
[A,~,c,dd] = LFR(n,d,mu,gamma, gamma_c, d_min);
flag = sum(sum(isnan(A))) > 0;
end
toc

%% DIFFUSIONE DI OPINIONE
epsilon = 0.12; % Livello di fiducia
% GENERAZIONE DELLA MATRICE DI TRUSTINESS
W = trustiness(A);

opin = rand(n,1); % OPINIONI INIZIALI
confidence = epsilon*ones(n,1); % LIVELLI DI CONFIDENZA
tic
[opinionHistory, it] = HK(A, W, opin, confidence); % DIFFUSIONE
toc


%% Grafici  consenso, polarizzazione e frammentazione dell'opinione
dinamicaOpinione(opinionHistory,0)