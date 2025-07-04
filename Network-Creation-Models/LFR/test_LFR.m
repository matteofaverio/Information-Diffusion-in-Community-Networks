clc,clear,close all
gamma = 2+1e-12;
d_min = 2;
n = 10000;
d = 10;
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
n = 1000;
d_min = 10;
d_max = 20;
gamma_c = 3;
[S,N] = powerLaw_communities(n,d_min,d_max,gamma_c);
figure;
histogram(S, 'Normalization','pdf', 'EdgeColor','none');
title('Istogramma dei campioni ~ Power Law troncata');
xlabel('x');
ylabel('Densità di probabilità stimata');

%% 
clc, clear, close all
n = 200;
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.9;
[A,AA,c,h,L,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
%%
G = graph(A,'upper','omitselfloops');
c = c(:);
comList = unique(c);          % etichette comunità
K = numel(comList);

XY = zeros(numnodes(G),2);    % coord finali di tutti i nodi
theta = linspace(0, 2*pi, K+1);   % K angoli equispaziati
R = 10;                          % raggio del cerchio (più alto → comunità più lontane)

for k = 1:K
    idx = (c == comList(k));      % nodi della k-esima comunità
    Gk  = subgraph(G, idx);       % sottografo

    % layout locale
    pk = plot(Gk,'Layout','force','Iterations',150,'UseGravity',false);
    x  = pk.XData;  y = pk.YData;

    % ridimensiona per tenerli in ~[-1,1]²
    s = 1 / max(range([x(:); y(:)])) * 1.6;
    x = x * s;  y = y * s;

    % offset sul cerchio
    off = [R*cos(theta(k)), R*sin(theta(k))];
    XY(idx,1) = x + off(1);
    XY(idx,2) = y + off(2);
end

% grafico finale
figure
h = plot(G, ...
    'XData',        XY(:,1), ...
    'YData',        XY(:,2), ...
    'MarkerSize',   5, ...
    'EdgeAlpha',    .15, ...
    'LineWidth',    .50, ...
    'NodeCData',    c);
colormap(lines(K));  
%%
% ---- INPUT -------------------------------------------------------------
G = graph(A,'upper','omitselfloops');
c = c(:);                       % vettore colonna
comList = unique(c);            % etichette comunità
K = numel(comList);

% ---- PARAMETRI CHE PUOI REGOLARE --------------------------------------
iterationsLocal = 300;          % più alto ⇒ layout interno più disteso
localScale      = 2.5;          % >1 ⇒ nodi interni più lontani
R               = 10;           % distanza fra i “centri” delle comunità
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
    'MarkerSize', 3, ...
    'EdgeAlpha',  0.25, ...
    'LineWidth',  .1);
colormap(lines(K));  
set(gca, 'XTick', [], 'YTick', []);

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
p = plot(graph(AA));
for i = 1:n
    highlight(p,i,'MarkerSize',log(c(i)+1),'NodeColor',color(c(i),:))
    p.NodeLabel = [];
end
figure(2)
p1 = plot(graph(A));
for i = 1:n
    highlight(p1,i,'MarkerSize',log(c(i)+1),'NodeColor',color(c(i),:))
    p1.NodeLabel = [];
end
