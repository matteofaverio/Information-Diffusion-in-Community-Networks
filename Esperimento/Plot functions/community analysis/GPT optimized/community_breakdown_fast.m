function results = community_breakdown_fast(folderName)
    % ----- PARAMETRI FISSI ---------------------------------------------
    epsilons = 0.01:0.01:0.25;     nEps   = numel(epsilons);
    nRuns    = 200;                % 200 matrici per file
    edges    = linspace(0,1,26);   % 25 bin  (usato dalla entropy)

    % ----- ALLOCAZIONE OUTPUT (numeric, NON cell) ----------------------
    Ent = zeros(nRuns, nEps);      % [run × eps]
    CI  = zeros(nRuns, nEps);
    BC  = zeros(nRuns, nEps);

    % ----- INFO SULLE COMUNITÀ (caricata UNA SOLA VOLTA) ---------------
    dataDir  = 'C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Data';
    muTok    = regexp(folderName,'(\d{3})','match','once');      % "050"…
    idxStr   = sprintf('%02d',(str2double(muTok)-50)/5+1);       % "01"…
    C        = load(fullfile(dataDir,"data_mu"+idxStr+".mat")).C; % 1×200 cell

    %   Pre-calcola indici interi per ogni nodo → comunità  (serve ad accumarray)
    commIdx  = cellfun(@(g) grp2idx(g), C, 'uni',0);
    nCommMax = max( cellfun(@max, commIdx) );

    % ----- ELENCO FILE DA ANALIZZARE -----------------------------------
    fileList = natsortfiles( dir(fullfile(folderName,'*.mat')) ); % ord. "naturale"
    assert(numel(fileList)==nEps,'Attesi %d file, trovati %d',nEps,numel(fileList))

    % ----- ELABORAZIONE PARALLELA --------------------------------------
    parfor e = 1:nEps
        S     = load( fullfile(folderName, fileList(e).name) );
        ops = struct2cell(S);
    ops = ops{1};           % 3×200 cell

        % vettori temporanei (parfor safe)
        Ent_e = zeros(nRuns,1);
        CI_e  = zeros(nRuns,1);
        BC_e  = zeros(nRuns,1);

        for r = 1:nRuns
            opin = double(ops{2,r})/1000;     % stato finale in [0,1]
            g    = commIdx{r};                % comunità indicizzate 1…k

            % -- STATISTICHE VETTORIALI ---------------------------------
            [ci, h, bc] = fast_Stats(opin, g, nCommMax, edges);

            % medie sulle comunità (se ti interessa la media)
            CI_e(r)  = mean(ci);
            Ent_e(r) = mean(h);
            BC_e(r)  = mean(bc);
        end

        CI (:,e) = CI_e;
        Ent(:,e) = Ent_e;
        BC (:,e) = BC_e;
        fprintf('ε = %.2f completato\n', epsilons(e))
    end

    % ----- STRUCT DI OUTPUT --------------------------------------------
    results.ShannonEntropy        = Ent;
    results.ConcentrationIndex    = CI;
    results.BimodalityCoefficient = BC;
end
