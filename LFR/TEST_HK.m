%% GENERAZIONE DEL NETWORK

clc; clear; close all;

% PARAMETRI INIZIALI
n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;

t = 1;
TEST = cell(t,t);
mu = 0.8; %linspace(0.3,0.999,t);
epsilon = 0.15; % linspace(0.05,0.2,t);

for i = 1:t

    % GENERAZIONE DELLA RETE
    flag = true;
    tic
    while flag
    [A,~,c,dd] = LFR2(n,d,mu(i),gamma, gamma_c, d_min);
    flag = sum(sum(isnan(A))) > 0;
    end
    toc
    
    % GENERAZIONE DELLA MATRICE DI TRUSTINESS
    W = trustiness(A);
    
    for j = 1:t

    % DIFFUSIONE DI OPINIONE
    opin = rand(n,1); % OPINIONI INIZIALI
    confidence = epsilon(j)*ones(n,1); % LIVELLI DI CONFIDENZA
    % agents = randi([1,n],0,1);
    tic
    [finalOpinions, it, opinionHistory] = HK(A, W, opin, confidence); % DIFFUSIONE
    toc

    TEST{i}{j} = opinionHistory;
    fprintf('mu = %d , eps = %d \n',mu(i),epsilon(j));

    end

    fprintf('\n');
end

%% MU FISSO, EPSILON VARIABILE

clc; clear; close all;

% PARAMETRI INIZIALI
n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;

n_var = 20;
n_tries = 10;
TEST = cell(n_var,n_tries);
mu = 0.75;
epsilon = linspace(0.01,0.3,t);


for i = 1:t

    [A,~,c,dd] = network_LFR(n,d,mu(i),gamma, gamma_c, d_min);
    
    W = trustiness(A);

    for j = 1:v
    
    opin = rand(n,1); 
    confidence = epsilon*ones(n,1); 
    [finalOpinions, it, opinionHistory] = HK(A, W, opin, confidence);
    TEST{i}{j} = opinionHistory;

    end
    fprintf('eps = %d \n\n ',epsilon(i));
end

%% EPSILON FISSO, MU VARIABILE

clc; clear; close all;

% PARAMETRI INIZIALI
n = 1000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;

n_var = 20;
n_tries = 10;
TEST = cell(n_var,n_tries);
mu = linspace(0.05,0.95,n_var);
epsilon = 0.18;

A = cell(n_var);
W = cell(n_var);

for i = 1:n_var
    [A{i},~,c,dd] = network_LFR(n,d,mu(i),gamma, gamma_c, d_min);
    W{i} = trustiness(A{i});
end

for i = 1:n_var
    for j = 1:n_tries
    
        opin = rand(n,1); 
        confidence = epsilon*ones(n,1); 
        [finalOpinions, it, TEST{i}{j}] = HK(A{i}, W{i}, opin, confidence);
    
    end
    fprintf('mu = %d \n\n',mu(i));
end



%% ANIMAZIONE VIDEO DELLA DIFFUSIONE

i = 1;
X = TEST{i}{j};

% Prepara la figura
figure;
set(gcf, 'Position', [100, 100, 600, 400]);  % Imposta dimensioni finestra (opzionale)

% (Opzionale) Se vuoi salvare il filmato come video, puoi preparare un oggetto VideoWriter:
writerObj = VideoWriter('convergenza.avi');
open(writerObj);

%X = opinionHistory;
T = size(X,2);
N = size(X,1);

% Inizializza un array di frame se vuoi usare la funzione movie in MATLAB
F(T) = struct('cdata',[],'colormap',[]);

x = 1:N;  % Asse x: gli indici dei valori

for t = 1:10:T
    % Prendi i valori alla colonna t
    valori = X(:,t);

    % Crea uno scatter dei punti, colorando ciascun punto in base al suo valore
    scatter(x, valori, 10, valori, 'filled');
    title(['Iterazione: ', num2str(t)]);
    xlabel('Indice del valore');
    ylabel('Valore');
    
    % Imposta i limiti dell’asse y in [0,1], se i tuoi valori sono sempre in quel range
    ylim([0, 1]);

    % Imposta la mappa di colori e la barra dei colori
    colormap(turbo);        % Scegli una mappa colori (es. 'jet', 'parula', 'hot', ecc.)
    colorbar;             % Mostra la barra di riferimento dei colori
    clim([0 1]);         % Se i valori sono in [0,1], fissa la scala dei colori

    %drawnow;              % Aggiorna istantaneamente il plot
    %pause(0.1);           % Metti in pausa per un breve intervallo (0.1 secondi)

    % Salva il frame per l'animazione o video
    F(t) = getframe(gcf);

    % Se vuoi salvare il video come .avi (o .mp4), potresti fare:
    writeVideo(writerObj, F(t));
end

% Se usi la funzione movie di MATLAB per rivedere l’animazione:
%movie(F, 1, 10);  % 1 ripetizione, 5 fps (fotogrammi al secondo)

% Se hai usato VideoWriter, alla fine chiudi l'oggetto:
close(writerObj);

%% HEATMAP DELLA DIFFUSIONE

n_tries = 1;

for i = 1:n_tries
    for j = 1:n_tries

    X = TEST{i}{j};
    [N, T] = size(X);
    
    % Numero di bin sui valori in [0,1]
    numBins = 50;
    edges = linspace(0, 1, numBins+1);
    
    % Calcola la matrice di densità H (numBins x T)
    H = zeros(numBins, T);
    for t = 1:T
        H(:, t) = histcounts(X(:, t), edges)';
    end
    
    % Crea la heatmap
    pos = (i - 1) * n_tries + j;
    subplot(n_tries,n_tries,pos)
    imagesc(H);
    axis xy;  % Per avere 0 in basso e 1 in alto
    colormap(spring);
    hcb = colorbar;
    hcb.Label.String = 'N° di nodi';
    title('Distribuzione delle opinioni durante le iterazioni');
    xlabel('Iterazioni');
    ylabel('Opinione');
    
    %--- PERSONALIZZAZIONE DELL'ASSE X (iterazioni) ---
    %Mostra tick ogni 20 iterazioni (o come preferisci)
    xticks = 0:50:T;  % ogni 20 iterazioni
    set(gca, 'XTick', xticks);
    set(gca, 'XTickLabel', string(xticks));
    
    %--- PERSONALIZZAZIONE DELL'ASSE Y (valori tra 0 e 1) ---
    %Calcola i centri dei bin
    binCenters = 0.5*(edges(1:end-1) + edges(2:end));
    
    %Scegli tick significativi per y, es: ogni 0.2
    yticks_vals = 0:0.2:1;
    %Trova gli indici più vicini nei binCenters
    [~, yticks] = min(abs(binCenters' - yticks_vals), [], 1);
    
    set(gca, 'YTick', yticks);
    set(gca, 'YTickLabel', string(yticks_vals));

    end
end
%%

% X: matrice (N x T), dove:
%   N = numero di valori (righe), T = numero di iterazioni (colonne)

[N, T] = size(X);

% Definiamo i bin per suddividere l'intervallo [0,1] in numBins intervalli
numBins = 20;
edges = linspace(0, 1, numBins + 1);                     % Estremi dei bin
binCenters = 0.5 * (edges(1:end-1) + edges(2:end));      % Centri dei bin

% Calcola la densità: H(i,j) = # di valori in bin i alla iterazione j
H = zeros(numBins, T);
for t = 1:T
    H(:, t) = histcounts(X(:, t), edges)';
end

% Prepara griglie per asse x (iterazioni) e y (valori)
[Xgrid, Ygrid] = meshgrid(1:T, binCenters);   % Xgrid = iterazioni, Ygrid = valori

% Grafico 3D della densità
figure;
surf(Xgrid, Ygrid, H, 'EdgeColor', 'none');   % 'EdgeColor' per una superficie liscia

% Etichette e stile
xlabel('Iterazione');
ylabel('Valore');
zlabel('Densità');
title('Distribuzione dei Valori nel Tempo (3D)');
colormap(turbo);             % Mappa colori per la superficie (facoltativo)
colorbar;                    % Colore aggiuntivo (utile, ma opzionale)
view(45, 30);                % Angolo di vista 3D
shading interp;              % Smussa la superficie
grid on;

%%

% Parametri
nEps   = n_var;  % numero di valori diversi di epsilon
nTests = n_tries;  % numero di test ripetuti per ogni epsilon
nBins  = 100; % numero di intervalli tra 0 e 1
myCell = TEST;

% Vettore dei bordi dei bin: 101 punti da 0 a 1
edges = linspace(0, 1, nBins+1);

% Preallocazione di una matrice che conterrà,
% per ogni epsilon (riga), l'istogramma normalizzato a 100 bin.
distMat = zeros(nEps, nBins);

% Calcolo delle distribuzioni
for i = 1:nEps
    % Raccolgo tutti i valori dei 10 test corrispondenti a epsilon i
    allVals = [];
    for j = 1:nTests
        M = myCell{i}{j};    % M è la matrice del test (dimensione n × c_j)
        allVals = [allVals; M(:)];  % Appiattisco e concateno
    end
    
    % Istogramma dei valori in 100 bin tra 0 e 1
    counts = histcounts(allVals, edges);
    
    % Normalizzo per avere una distribuzione di probabilità
    distMat(i,:) = counts / sum(counts);
end

% Supponendo di avere i valori effettivi di epsilon in un vettore 'epsVals'
% lungo 40. Se non lo hai, puoi usare semplicemente 1:40 come riferimento.
epsVals = mu;  % O il vero vettore dei 40 epsilon

% Per l'asse dei bin, prendiamo i centri di ogni bin (invece di usare i bordi)
binCenters = edges(1:end-1) + diff(edges)/2;

% Creiamo la griglia per il plot 3D.
%  - Asse X: i 40 valori di epsilon
%  - Asse Y: i 100 intervalli
%  - Asse Z: il valore della distribuzione
[X, Y] = meshgrid(epsVals, binCenters);
Z = distMat';  

% Soglia sotto cui considerare la distribuzione "bassa"
% threshold = 0.0001;
% Copia la Z originale
% Zplot = Z;
% Imposta a NaN i valori molto bassi (trasparente/bianco nel plot)
% Zplot(Zplot < threshold) = NaN;

% Plot in 3D
figure;
surf(X, Y, Z);
nColors = 256;
customMap = [ ...
    1.00 1.00 0.00;   % giallo
    1.00 0.00 0.00;   % rosso
    0.60 0.00 0.80;   % viola
    0.00 0.00 1.00];  % blu
cmap = interp1(linspace(0,1,size(customMap,1)), customMap, linspace(0,1,nColors));
colormap(cmap);
colorbar;  % facoltativo: aggiunge barra laterale dei colori
