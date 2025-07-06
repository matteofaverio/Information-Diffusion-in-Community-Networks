
function [ci, H, bc] = stats_corrette(op, grp, edges)
% Calcola CI, Entropia e BC in modo vettoriale e con la logica corretta.
% 'grp' deve essere un vettore di indici interi positivi (da grp2idx).
    
    % Dati validi (ignora eventuali NaN)
    valid_idx = ~isnan(op);
    op = op(valid_idx);
    grp = grp(valid_idx);

    nCommMax = max(grp);
    if isempty(nCommMax), nCommMax = 0; end % Gestisce il caso di input vuoto

    % 1) Concentration Index (delta=0.05)
    mu_ci  = accumarray(grp, op, [nCommMax 1], @mean, NaN);
    close_nodes = abs(op - mu_ci(grp)) <= 0.05;
    ci = accumarray(grp, close_nodes, [nCommMax 1], @mean, NaN);
    
    % 2) Shannon Entropy (CORRETTA)
    nBin   = numel(edges) - 1;
    nNode  = accumarray(grp, 1, [nCommMax 1]);
    
    op(op == 1) = 1 - 1e-9; % Evita caso limite di discretize
    bin    = discretize(op, edges);
    
    % 'accumarray' qui è più sicuro di hc'./nNode se qualche nNode è 0
    counts = accumarray([grp, bin], 1, [nCommMax, nBin]);
    p = counts ./ nNode; % nNode è già una colonna, broadcasting corretto
    
    % Logica corretta per evitare log(0)
    log_p = log(p);
    log_p(p == 0) = 0;
    
    % La somma è sulla seconda dimensione (sui bin)
    H_unnormalized = -sum(p .* log_p, 2);
    H = H_unnormalized / log(nBin);
    
    % 3) Bimodality Coefficient (con bias-correction, IDENTICO ALL'ORIGINALE)
    % Usiamo splitapply con le funzioni native di MATLAB, è veloce e garantisce
    % la correttezza del calcolo (inclusa la bias correction).
    % fcn_bc = @(x) (kurtosis(x, 1) / (skewness(x, 1)^2 + 1));
    fcn_bc = @(x) ( (kurtosis(x, 1) - 2) ./ (skewness(x, 1).^2 + 1) );
    bc = splitapply(fcn_bc, op, grp);
    
    % Assicura che l'output sia un vettore colonna della dimensione giusta,
    % riempiendo con NaN i gruppi non presenti in 'grp' ma inferiori a nCommMax.
    temp_bc = NaN(nCommMax, 1);
    unique_groups = unique(grp);
    temp_bc(unique_groups) = bc;
    bc = temp_bc;
end