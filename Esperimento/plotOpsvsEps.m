function [] = plotOpsvsEps(TEST,n_var,n_tries,)
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
end

