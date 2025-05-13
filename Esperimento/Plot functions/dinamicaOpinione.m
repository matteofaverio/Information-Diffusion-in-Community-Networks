function [] = dinamicaOpinione(X,flag)
% X: matrice (N x T), dove:
%   N = numero di valori (righe), T = numero di iterazioni (colonne)
if nargin == 1
    flag = 0;
end
[N, T] = size(X);

% Definiamo i bin per suddividere l'intervallo [0,1] in numBins intervalli
numBins = 50;
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
opdin = surf(Xgrid, Ygrid, H);   % 'EdgeColor' per una superficie liscia
set(opdin, 'EdgeColor', 'none')
% Etichette e stile
xlabel('Iterazione');
ylabel('Valore Opinione');
zlabel('Frequenza');
title('Distribuzione dei Valori nel Tempo (3D)');
colormap(turbo);             % Mappa colori per la superficie (facoltativo)
colorbar;                    % Colore aggiuntivo (utile, ma opzionale)
view(45, 30);                % Angolo di vista 3D
shading interp;              % Smussa la superficie
grid on;

if flag == 1
    n = convergence(X,1e-04);
    hold on 
    x0 = n;            
    yl = ylim;          
    zl = zlim;         
    [yp, zp] = meshgrid(linspace(yl(1), yl(2), 2), ...
                        linspace(zl(1), zl(2), 2));
    xp = x0 * ones(size(yp));
    hPlane = surf(xp, yp, zp);
    set(hPlane, 'FaceAlpha', 0.8, 'EdgeColor', 'w');
end
end

