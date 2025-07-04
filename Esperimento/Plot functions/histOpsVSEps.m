function histOpsVSEps(folderName)
% plotEpsilonDistribution  Istogramma 3D dello stato finale vs epsilon
%   plotEpsilonDistribution(folderName) carica tutti i .mat in folderName,
%   estrae dalle 200 matrici in ciascun file l'ultima colonna (stato finale),
%   calcola per ogni ε la frequenza normalizzata sui bin in [0,1], media su
%   tutti i test, e infine la visualizza con parallelepipedi 3D (bar3).
%   Mantiene la scala colori per l'asse z.

    % --- Parametri ---
    epsilons   = 0.01:0.01:0.25;      % 25 valori di epsilon
    nEps       = numel(epsilons);
    nBins      = 25;                  % numero di bin per l'istogramma
    binEdges   = linspace(0,1,nBins+1);
    binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;

    % --- Preallocazione ---
    histDist = zeros(nBins, nEps);

    % --- Elenco e ordinamento file ---
    files = dir(fullfile(folderName,'*.mat'));
    %files = files(~strncmp({files.name}, '._', 2));
    if numel(files)~=nEps
        error('Trovati %d file, ma ci aspettavamo %d (numero di epsilon)', ...
              numel(files), nEps);
    end
    names = sort({files.name});

    % --- Ciclo su ciascun file / valore di epsilon ---
    for i = 1:nEps
        
        S = load(fullfile(folderName, names{i}));
        fn = fieldnames(S);
        C  = S.(fn{1}); % cell array 3x200

        sumDist = zeros(nBins,1);
        for j = 1:200
            col  = double(C{2,j})/1000;   % scala
            d   = histcounts(col, binEdges, 'Normalization','probability');
            sumDist = sumDist + d';
        end
        histDist(:,i) = sumDist ./ numel(C);
    end

    % --- Colormap personalizzata ---
    nColors = 256;
    customMap = [ ...
        1.00 1.00 0.00;   % giallo
        1.00 0.00 0.00;   % rosso
        0.60 0.00 0.80;   % viola
        0.00 0.00 1.00];  % blu
    cmap = interp1(linspace(0,1,size(customMap,1)), customMap, linspace(0,1,nColors));

    % --- Plot 3D con parallelepipedi (bar3), con assi scambiati ---
    figure;
    hBars = bar3(histDist', 'detached');  % histDist' è nEps x nBins
    hold on;
    for k = 1:length(hBars)
        zdata = hBars(k).ZData;
        hBars(k).CData = zdata;
        hBars(k).FaceColor = 'interp';
        hBars(k).EdgeColor = 'k';
        hBars(k).LineWidth = 0.5;
    end

    % etichette e vista
    xlabel('Media opinioni finali','FontSize',12);
    ylabel('\epsilon','FontSize',12);
    zlabel('Frequenza media','FontSize',12);
    title('Distribuzione dello stato finale vs \epsilon (Istogramma 3D)','FontSize',14);
    view(230,30);
    grid on;

    % scala colori
    colormap(cmap);
    clim([0, max(histDist(:))]);

        % … tutto il resto del tuo codice …

    % regola ticks: 
    %  - X = opinioni [0:0.1:1]
    %  - Y = ε selezionati
    valoriOpinione = 0:0.25:1;                       % 0,0.1,…,1.0
    % mappa ciascun valore [0,1] sulla posizione di bar3 (coordinate 0.5–nBins+0.5)
    posX = valoriOpinione * nBins + 0.5;           
    etichetteX = arrayfun(@(v)sprintf('%.1f',v), valoriOpinione, 'UniformOutput',false);

    % scegli su quali indici di ε vuoi i tick (esempio 1,5,10,15,20,25)
    idxY = [1,5,10,15,20,25];
    posY = idxY;                                    % bar3 posiziona i gruppi in y=1:nEps
    etichetteY = arrayfun(@(i)sprintf('%.2f',epsilons(i)), idxY, 'UniformOutput',false);
    ylim([0.5, nEps+0.5]);

    % applica
    set(gca, ...
    'XLim',   [0.5,        nBins+0.5], ...
    'YLim',   [0.5,        nEps+0.5], ...
    'ZLim',   [0,        0.4], ...
    'XTick',  posX, ...
    'XTickLabel', etichetteX, ...
    'YTick',  posY, ...
    'YTickLabel', etichetteY, ...
    'Box',    'on', ...
    'ZGrid',  'off', ...
    'Layer',  'top');


end
