%% Generazione Rete
clc, clear, close all
n = 10000;      
gamma = 3;
gamma_c = 3;
d = 55;
d_min = 30;

t = 1;
TEST = cell(t,t);
mu = 0.99; %linspace(0.3,0.999,t);

flag = true;
tic
while flag
[A,~,c,dd] = LFR(n,d,mu,gamma, gamma_c, d_min);
flag = sum(sum(isnan(A))) > 0;
end
toc

%% DIFFUSIONE DI OPINIONE
epsilon = 0.13; % Livello di fiducia
% GENERAZIONE DELLA MATRICE DI TRUSTINESS
W = trustiness(A);

opin = rand(n,1); % OPINIONI INIZIALI
confidence = epsilon*ones(n,1); % LIVELLI DI CONFIDENZA
tic
[opinionHistory, it] = HK(A, W, opin, confidence); % DIFFUSIONE
toc
%% Grafici  consenso, polarizzazione e frammentazione dell'opinione
dinamicaOpinione(opinionHistory,0)
%%
% ---- INPUT -------------------------------------------------------------
G = graph(A,'upper','omitselfloops');
c = c(:);                       % vettore colonna
comList = unique(c);            % etichette comunità
K = numel(comList);

% ---- PARAMETRI CHE PUOI REGOLARE --------------------------------------
iterationsLocal = 300;          % più alto ⇒ layout interno più disteso
localScale      = 1;          % >1 ⇒ nodi interni più lontani
R = 20;          % distanza fra i "centri" delle comunità
% -----------------------------------------------------------------------

XY    = zeros(numnodes(G),2);   % coord. globali
theta = linspace(0,2*pi,K+1);   % angoli sul cerchio

for k = 1:K
    idx = (c == comList(k));    % nodi della k-esima comunità
    Gk  = subgraph(G,idx);

    % 1) --- ottieni un layout FORCE *solo* per la comunità --------------
    f   = figure('Visible','off');             % evita aprire finestre
    pk  = plot(Gk, ...
               'Layout','force', ...
               'Iterations',iterationsLocal);  % ← qui spingi la repulsione
    x = pk.XData;  y = pk.YData;
    close(f);                                  % chiudi la figura invisibile

    % 2) --- scala per “esplodere” la comunità ---------------------------
    s = localScale / max(range([x; y]));
    x = x * s;
    y = y * s;

    % 3) --- offset sul cerchio (o dove preferisci) ----------------------
    off = [R*cos(theta(k)), R*sin(theta(k))];
    XY(idx,1) = x + off(1);
    XY(idx,2) = y + off(2);
end

% ---- GRAFICO FINALE ----------------------------------------------------
figure
h = plot(G, ...
    'XData',      XY(:,1), ...
    'YData',      XY(:,2), ...
    'NodeCData',  c, ...
    'MarkerSize', 2, ...
    'EdgeAlpha',  .15, ...
    'LineWidth',  .15);
colormap(lines(K));  
ax = gca;
ax.XTick = [];
ax.YTick = [];
title("Visualizzazione delle comunita'",'Interpreter','latex')

%% Visualizzazione consenso nelle comunità
% --- INPUT ---------------------------------------------------------------
% A       : matrice di adiacenza (NxN)
% c       : vettore (Nx1) con comunità di ogni nodo
% opinion : vettore (Nx1) con valore di opinione in [0,1] per ogni nodo
% -------------------------------------------------------------------------

opinion = opinionHistory(:,end);
G        = graph(A,'upper','omitselfloops');
c        = c(:);
comList  = unique(c);
K        = numel(comList);

% Parametri di layout
iterationsLocal = 300;          % più alto ⇒ layout interno più disteso
localScale      = 2;          % >1 ⇒ nodi interni più lontani
R = 20;          % distanza fra i "centri" delle comunità

% Preallocazione
XY    = zeros(numnodes(G),2);
theta = linspace(0,2*pi,K+1);

% --- Calcolo delle coordinate (Strategia B modificata) ------------------
for k = 1:K
    idx = (c == comList(k));
    Gk  = subgraph(G, idx);

    % 1) layout “force” solo per la comunità
    f  = figure('Visible','off');
    pk = plot(Gk, 'Layout','force', 'Iterations',iterationsLocal);
    x  = pk.XData;  y = pk.YData;
    close(f);

    % 2) scala interna
    s = localScale / max(range([x; y]));
    x = x * s;
    y = y * s;

    % 3) offset sul cerchio
    off = [R*cos(theta(k)), R*sin(theta(k))];
    XY(idx,1) = x + off(1);
    XY(idx,2) = y + off(2);
end

% --- Grafico finale ----------------------------------------------------
figure
h = plot(G, ...
    'XData',      XY(:,1), ...
    'YData',      XY(:,2), ...
    'NodeCData',  opinion, ...      % <-- qui colore secondo l'opinione
    'MarkerSize', 3, ...
    'EdgeAlpha',  .1, ...
    'LineWidth',  .25);

% togli tick
set(gca, 'XTick', [], 'YTick', []);
axis('square')

% colormap continua e scala [0,1]
colormap("cool");
clim([0 1]);

% colorbar con etichette da 0 a 1
hcb = colorbar;
hcb.Ticks      = [0 0.5 1];
hcb.TickLabels = {'0','0.5','1'};
ylabel(hcb, 'Opinione', 'Interpreter','latex','FontSize',15);

title("Visualizzazione del consenso nelle comunita'",'Interpreter','latex','fontSize',15)
%% Grafici evoluzione opinione
main_folder = 'Results';
iterazioni(main_folder)
