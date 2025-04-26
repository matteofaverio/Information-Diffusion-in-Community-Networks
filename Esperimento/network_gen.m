%% ------------------------------------------------------------------------
% Aggiunge "Network", "Data" e "Results" al search-path indipendentemente
% da come avvii main.m (Run/F5, riga per riga, da Command Window, ecc.).
% ------------------------------------------------------------------------

% 1) Trova il percorso reale di questo script
if usejava('desktop')                          % sessione desktop
    doc = matlab.desktop.editor.getActive;     % documento aperto in Editor
    if ~isempty(doc) && ~isempty(doc.Filename) % file salvato su disco
        scriptPath = doc.Filename;             % ← percorso reale
    end
end

% 2) Se non arriva dall'Editor, prova con mfilename (funziona in CLI/run)
if ~exist('scriptPath','var') || isempty(scriptPath)
    tmp = mfilename('fullpath');
    if ~isempty(tmp)
        scriptPath = tmp;
    else
        % 3) Ultima risorsa: assumi che main.m sia nella Current Folder
        scriptPath = fullfile(pwd,'main.m');
    end
end

projDir  = fileparts(scriptPath);              % cartella del progetto
subDirs  = {'Network','Data','Results'};       % quelle da aggiungere

for k = 1:numel(subDirs)
    thisDir = fullfile(projDir, subDirs{k});
    if exist(thisDir,'dir')
        addpath(thisDir);                      % oppure addpath(genpath(thisDir))
    else
        warning('La cartella "%s" non esiste.', thisDir);
    end
end

clear doc tmp scriptPath projDir subDirs k     % pulizia variabili

%% GENERATION

clear
clc

% Parametri
n       = 10000;
gamma   = 3;
gamma_c = 2;
d       = 10;
d_min   = 7;

tests   = 200; 
mu      =  [0.51 0.56 0.62 0.64 0.70 0.76 0.78 0.84 0.88 0.92];

% Cartella di output
outDir = 'Output_gen';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for i = 1:length(mu)

    % Nome del file per questo valore di mu
    outFile = fullfile(outDir, sprintf('data_mu%02d.mat', i));
    if exist(outFile, 'file')
        delete(outFile);
    end

    % 1) Preallocazione su disco
    M       = cell(1, tests);   % per le matrici A
    C       = cell(1, tests);   % per i vettori c
    AVG_mu = zeros(1, tests);  % per i risultati AVG_MU(i,:)
    save(outFile, 'M', 'C', 'AVG_mu', '-v7.3');
    clear M C AVG_mu;

    % 2) Apre il MAT in scrittura parziale
    matObj = matfile(outFile, 'Writable', true);

    % 3) Loop interno: genera A e c, calcola, salva e libera RAM
    for j = 1:tests
        [A, ~, c, ~] = LFR2(n, d, mu(i), gamma, gamma_c, d_min);
      
        % calcolo AVG_MU(i,j)
        Mcomm           = (c == c');
        sameCommCounts  = sum(A .* Mcomm, 2);
        degs            = sum(A, 2);
        fracs           = sameCommCounts ./ degs;
        avg_val         = mean(fracs);

        % riduzione della matrice a una sparsa
        A = sparse(A);

        % scrive parzialmente su disco
        matObj.M(1, j)       = {A};
        matObj.C(1, j)       = {c};
        matObj.AVG_mu(1, j) = avg_val;

        clear A c
        fprintf('  [%2d,%3d] salvato in %s\n', i, j, outFile);
    end

    fprintf('Completato file %s (mu=%.2f)\n\n', outFile, mu(i));
end

%% REWIRING

mu_target = linspace(0.5,0.95,10);

for i = 1%:length(mu_values)

    % Costruisco il nome del file
    fname = fullfile('Output_gen', sprintf('data_mu%02d.mat', i));
    if ~exist(fname, 'file')
        warning('File non trovato: %s', fname);
        continue;
    end

    % Apro in scrittura parziale
    matObj = matfile(fname, 'Writable', true);

    % Leggo tutto AVG_mu in memoria
    avg_mu = matObj.AVG_mu;  

    % Trovo gli indici da sistemare
    idx_fix = find(avg_mu > mu_target(i)+0.005);

    if isempty(idx_fix)
        fprintf('Nessuna matrice da risistemare in %s\n', fname);
        continue;
    end

    fprintf('Riparazione in %s: %d matrici su 200\n', fname, numel(idx_fix));
    for j = idx_fix
        % Carico la matrice e il vettore comunità
        tmpM = matObj.M(1, j);
        A    = tmpM{1};
        tmpC = matObj.C(1, j);
        c    = tmpC{1};
        
        % Trasformo la matrice sparsa e applico il rewiring
        A = full(A);
        A_new = rewiring(A, c, mu_target(i),1e4);

        % Riscrivo la matrice corretta
        matObj.M(1, j) = {A_new};

        fprintf('  → Matrice %d corretta\n', j);
    end
end
