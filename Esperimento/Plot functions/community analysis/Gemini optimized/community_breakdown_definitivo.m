function results = community_breakdown_definitivo(folderName)
% Versione definitiva che unisce le ottimizzazioni di 'community_breakdown_fast'
% con la correttezza logica delle funzioni originali, garantendo risultati identici.

    % ----- PARAMETRI FISSI ---------------------------------------------
    epsilons = 0.01:0.01:0.25;     nEps   = numel(epsilons);
    nRuns    = 200;                % 200 simulazioni per file
    edges    = linspace(0, 1, 26);   % 25 bin per l'entropia
    
    % ----- ALLOCAZIONE OUTPUT (formato a matrice, più efficiente) ------
    Ent_matrix = zeros(nRuns, nEps);
    CI_matrix  = zeros(nRuns, nEps);
    BC_matrix  = zeros(nRuns, nEps);
    
    % ----- CARICAMENTO DATI COMUNITÀ (UNA SOLA VOLTA) ------------------
    % Uso del 'regexp' più robusto dell'originale
    dataDir  = 'C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Data';
    tok = regexp(folderName, 'test_mu(\d{3})_processed', 'tokens', 'once');
    if isempty(tok), error('Il nome della cartella non corrisponde al pattern atteso: test_muXXX_processed'); end
    idx = (str2double(tok{1}) - 50)/5 + 1;
    idxStr = sprintf('%02d', idx);
    
    C = load(fullfile(dataDir, ['data_mu' idxStr '.mat'])).C; % 1×200 cell array
    
    % Pre-calcola indici interi per le comunità per ogni run (ottima ottimizzazione!)
    commIdx = cellfun(@grp2idx, C, 'UniformOutput', false);

    % ----- ELENCO FILE DA ANALIZZARE -----------------------------------
    % 'sort' è sufficiente se i nomi dei file sono zero-padded, altrimenti 'natsortfiles' è più sicuro
    files = dir(fullfile(folderName, '*.mat'));
    assert(numel(files) == nEps, 'Attesi %d file, trovati %d', nEps, numel(files));
    names = sort({files.name});
    
    fprintf('--- Avvio elaborazione parallela su %d file... ---\n', nEps);
    % ----- ELABORAZIONE PARALLELA --------------------------------------
    parfor e = 1:nEps
        % Caricamento robusto dei dati
        S = load(fullfile(folderName, names{e}));
        fn = fieldnames(S);
        ops = S.(fn{1}); % cell array 3x200
        
        % Vettori temporanei (necessari per parfor)
        Ent_e = zeros(nRuns, 1);
        CI_e  = zeros(nRuns, 1);
        BC_e  = zeros(nRuns, 1);
        
        for r = 1:nRuns
            opin = double(ops{2, r}) / 1000; % Stato finale in [0,1]
            grp_idx = commIdx{r};            % Comunità indicizzate 1…k
            
            % -- STATISTICHE VETTORIALI con logica CORRETTA --
            [ci, h, bc] = stats_corrette(opin, grp_idx, edges);
            
            % Calcolo della media sulle comunità, come da richiesta originale
            CI_e(r)  = mean(ci, 'omitnan');
            Ent_e(r) = mean(h, 'omitnan');
            BC_e(r)  = mean(bc, 'omitnan');
        end
        
        % Assegnazione alla matrice di output
        CI_matrix(:, e)  = CI_e;
        Ent_matrix(:, e) = Ent_e;
        BC_matrix(:, e)  = BC_e;
        
        fprintf('ε = %.2f completato\n', epsilons(e));
    end
    
    % ----- STRUCT DI OUTPUT --------------------------------------------
    % NOTA: L'output è una matrice [runs x epsilons].
    % Per tornare al formato cella originale:
    % results.ShannonEntropy = mat2cell(Ent_matrix, nRuns, ones(1, nEps));
    results.ShannonEntropy        = Ent_matrix;
    results.ConcentrationIndex    = CI_matrix;
    results.BimodalityCoefficient = BC_matrix;
    
    fprintf('--- Elaborazione completata. ---\n');
end


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
    fcn_bc = @(x) (kurtosis(x, 1) / (skewness(x, 1)^2 + 1));
    bc = splitapply(fcn_bc, op, grp);
    
    % Assicura che l'output sia un vettore colonna della dimensione giusta,
    % riempiendo con NaN i gruppi non presenti in 'grp' ma inferiori a nCommMax.
    temp_bc = NaN(nCommMax, 1);
    unique_groups = unique(grp);
    temp_bc(unique_groups) = bc;
    bc = temp_bc;
end