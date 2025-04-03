%% GENERAZIONE DEL NETWORK

clc; clear; close all;

% PARAMETRI INIZIALI
n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;
mu = 0.8;

% GENERAZIONE DELLA RETE
tic
[A,~,c,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min);
toc

%% GENERAZIONE DELLA MATRICE DI TRUSTINESS
tic
W = trustiness(A);
toc

%% DIFFUSIONE DI OPINIONE

opin = rand(n,1); % OPINIONI INIZIALI
confidence = 0.1*ones(n,1); % LIVELLI DI CONFIDENZA
tic
[finalOpinions, it, opinionHistory] = HK(A, W, opin, confidence); % DIFFUSIONE
toc


%% ANIMAZIONE VIDEO DELLA DIFFUSIONE

% X è una matrice N x T
% Ad esempio, per test potresti generarla così:
% N = 10; T = 50;
% X = rand(N, T);   % Valori iniziali uniformemente distribuiti in [0,1]
% (Nel tuo caso X è già definita e si aggiorna iterativamente.)

% Prepara la figura
figure;
set(gcf, 'Position', [100, 100, 600, 400]);  % Imposta dimensioni finestra (opzionale)

% (Opzionale) Se vuoi salvare il filmato come video, puoi preparare un oggetto VideoWriter:
%writerObj = VideoWriter('convergenza.avi');
%open(writerObj);

X = opinionHistory;
T = size(X,2);
N = size(X,1);

% Inizializza un array di frame se vuoi usare la funzione movie in MATLAB
F(T) = struct('cdata',[],'colormap',[]);

x = 1:N;  % Asse x: gli indici dei valori

for t = 1:5:T
    % Prendi i valori alla colonna t
    valori = X(:,t);

    % Crea uno scatter dei punti, colorando ciascun punto in base al suo valore
    scatter(x, valori, 15, valori, 'filled');
    title(['Iterazione: ', num2str(t)]);
    xlabel('Indice del valore');
    ylabel('Valore');
    
    % Imposta i limiti dell’asse y in [0,1], se i tuoi valori sono sempre in quel range
    ylim([0, 1]);

    % Imposta la mappa di colori e la barra dei colori
    colormap(jet);        % Scegli una mappa colori (es. 'jet', 'parula', 'hot', ecc.)
    colorbar;             % Mostra la barra di riferimento dei colori
    clim([0 1]);         % Se i valori sono in [0,1], fissa la scala dei colori

    drawnow;              % Aggiorna istantaneamente il plot
    pause(0.1);           % Metti in pausa per un breve intervallo (0.1 secondi)

    % Salva il frame per l'animazione o video
    F(t) = getframe(gcf);

    % Se vuoi salvare il video come .avi (o .mp4), potresti fare:
    %writeVideo(writerObj, F(t));
end

% Se usi la funzione movie di MATLAB per rivedere l’animazione:
%movie(F, 1, 10);  % 1 ripetizione, 5 fps (fotogrammi al secondo)

% Se hai usato VideoWriter, alla fine chiudi l'oggetto:
close(writerObj);

%% SCATTER PLOT DELLA DIFFUSIONE

figure;
hold on;
X = opinionHistory;
T = size(X,2);
N = size(X,1);

% Disegna tutti i punti in un singolo grafico
for t = 1:T
    % Estrai i valori alla colonna t
    yvals = X(:, t);
    % La coordinata x è fissa a t per tutti i valori di questa colonna
    xvals = t * ones(N, 1);
    % Scatter colorato in base al valore effettivo (yvals)
    scatter(xvals, yvals, 5, yvals, 'filled'); 
end

% Impostazioni dell'asse e della mappa colori
xlabel('Iterazione');
ylabel('Valore');
ylim([0, 1]);          % Se sai che i valori sono in [0,1]
clim([0, 1]);         % Fissa la scala colori da 0 a 1
colormap(jet);         % Scegli la mappa (es. 'jet', 'hot', ecc.)
colorbar;              % Barra dei colori
title('Convergenza dei Valori');

hold off;

%% HEATMAP DELLA DIFFUSIONE

% Supponiamo che X sia N x T
[N, T] = size(X);

% Numero di bin sui valori in [0,1]
numBins = 25;
edges = linspace(0, 1, numBins+1);

% Calcola la matrice di densità H (numBins x T)
H = zeros(numBins, T);
for t = 1:T
    H(:, t) = histcounts(X(:, t), edges)';
end

% Crea la heatmap
figure;
imagesc(H);
axis xy;  % Per avere 0 in basso e 1 in alto
colormap(hot);
hcb = colorbar;
hcb.Label.String = 'N° di nodi';
title('Distribuzione delle opinioni durante le iterazioni');
xlabel('Iterazioni');
ylabel('Opinione');

% --- PERSONALIZZAZIONE DELL'ASSE X (iterazioni) ---
% Mostra tick ogni 20 iterazioni (o come preferisci)
xticks = 0:50:T;  % ogni 20 iterazioni
set(gca, 'XTick', xticks);
set(gca, 'XTickLabel', string(xticks));

% --- PERSONALIZZAZIONE DELL'ASSE Y (valori tra 0 e 1) ---
% Calcola i centri dei bin
binCenters = 0.5*(edges(1:end-1) + edges(2:end));

% Scegli tick significativi per y, es: ogni 0.2
yticks_vals = 0:0.2:1;
% Trova gli indici più vicini nei binCenters
[~, yticks] = min(abs(binCenters' - yticks_vals), [], 1);

set(gca, 'YTick', yticks);
set(gca, 'YTickLabel', string(yticks_vals));
