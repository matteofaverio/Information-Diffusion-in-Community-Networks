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
addpath(genpath('C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks'));


n_tests = 100;
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
nmis = zeros(n_tests,mm,gg,bb,dd);

for i = 1:n_tests
    for m = 1:mm
        for g = 1:gg
            for b = 1:bb
                for d = 1:dd
                    [A(:,:,m,g,b,d),AA,c(:,m,g,b,d),degrees] = network_LFR(n,d_mean(d),mu(m),gamma(g), beta(b), d_min(d));
                    if sum(sum(isnan(A))) > 0
                        nmis(i,m,g,b,d) = nan;
                        continue
                    end
                    Q = community_louvain(A(:,:,m,g,b,d));
                    nmis(i,m,g,b,d) = nmi(c(:,m,g,b,d),Q);
                end
            end
        end
    end
    fprintf('iterazione %d\n',i)
end

save('nmiresults.mat','nmis');
%%



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

%% 📝 MATLAB Cheat Sheet – Animare un plot

% Questo script mostra diversi metodi per animare grafici in MATLAB.
% Ogni sezione è indipendente. Puoi eseguire tutto o solo una parte.

%% 1️⃣ Animare un punto su una sinusoide

x = linspace(0, 2*pi, 100);
y = sin(x);

figure;
h = plot(x(1), y(1), 'ro', 'MarkerSize', 10); % punto iniziale (rosso)
xlim([0 2*pi]);
ylim([-1.2 1.2]);

for k = 2:length(x)
    set(h, 'XData', x(k), 'YData', y(k)); % aggiorna posizione del punto
    pause(0.05); % controlla la velocità dell'animazione
end

%% 2️⃣ Disegnare una curva progressivamente (linea che si estende)

x = linspace(0, 2*pi, 100);
y = sin(x);

figure;
h = plot(NaN, NaN, 'b-', 'LineWidth', 2); % linea vuota all'inizio
xlim([0 2*pi]);
ylim([-1.2 1.2]);

for k = 1:length(x)
    set(h, 'XData', x(1:k), 'YData', y(1:k)); % aggiorna la curva fino al punto k
    pause(0.03); % tempo tra i frame
end

%% 3️⃣ Evitare flickering e migliorare fluidità

% drawnow forza l’aggiornamento del grafico immediatamente
% utile quando si fanno più aggiornamenti rapidi

x = linspace(0, 2*pi, 100);
y = cos(x);

figure;
h = plot(NaN, NaN, 'g-', 'LineWidth', 2);
xlim([0 2*pi]);
ylim([-1.2 1.2]);

for k = 1:length(x)
    set(h, 'XData', x(1:k), 'YData', y(1:k));
    drawnow; % forza refresh del grafico
    pause(0.02);
end

%% 4️⃣ Animare un punto che segue una traiettoria (es. cerchio)

t = 0:0.1:10;
x = cos(t);
y = sin(t);

figure;
plot(x, y, 'k--'); % traiettoria di riferimento (opzionale)
hold on;
punto = plot(NaN, NaN, 'ro', 'MarkerSize', 10); % punto mobile

for k = 1:length(t)
    set(punto, 'XData', x(k), 'YData', y(k));
    drawnow;
    pause(0.05);
end

%% ℹ️ Note utili

% - set(obj, 'Prop', value): aggiorna proprietà di un oggetto grafico
% - pause(t): ferma l’esecuzione per t secondi (es. pause(0.05))
% - drawnow: forza aggiornamento immediato del plot (più fluido)
% - hold on/off: permette di disegnare più oggetti sullo stesso grafico
% - xlim / ylim: imposta limiti degli assi, utile per evitare flickering

% ✅ Fine cheat sheet
