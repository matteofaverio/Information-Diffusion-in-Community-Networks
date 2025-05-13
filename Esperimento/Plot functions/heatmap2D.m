function [] = heatmap2D(X)
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
    imagesc(H);
    axis xy;  % Per avere 0 in basso e 1 in alto
    colormap(summer);
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

