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

%% TEST


N = 1e4;
mu = linspace(0.5,0.95,10);
tests = 100; 
epsilons = linspace(0.01,0.25,25)';

for i = 10
    fname = fullfile('output_rewired', sprintf('data_mu%02d.mat', i));
    if ~exist(fname, 'file')
        warning('File not found: %s', fname);
        continue;
    end
    S = load(fname, 'M');

    folderPath = fullfile('output_rewired', sprintf('TEST_mu%g', mu(i)));
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    for j = 1:25
        
        opinions_j = cell(1, tests);
        it_vec = zeros(1, tests);
        tic
        for k = 1:tests
            A = full(S.M{k+100});
            W = trustiness(A);
            opin0 = rand(N, 1);
            [opinHistory, it] = HK(A, W, opin0, epsilons(j));
            scale = 1000;
            opinions_j{k} = uint16(round(opinHistory * scale));
            it_vec(k) = it;
        end
        toc

        avg_it = round(mean(it_vec));
        mu_str  = sprintf('%02d', round(mu(i) * 10));
        eps_str = sprintf('%03d', round(epsilons(j) * 100));
        varName = sprintf('opinions_mu%s_eps%s', mu_str, eps_str);

        save( fullfile(folderPath, [varName '.mat']), 'opinions_j', '-v7.3' );
        fprintf('μ = %g | ε = %g | it = %d\n\n', mu(i), epsilons(j), avg_it);
    end
end
   