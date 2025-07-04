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
subDirs  = {'Network&Test functions','Data','Results','Plot functions'};       % quelle da aggiungere

for k = 1:numel(subDirs)
    thisDir = fullfile(projDir, subDirs{k});
    if exist(thisDir,'dir')
        addpath(thisDir);                      % oppure addpath(genpath(thisDir))
    else
        warning('La cartella "%s" non esiste.', thisDir);
    end
end

clear doc tmp scriptPath projDir subDirs k     % pulizia variabili

%%  primo test da cancellare
n = 2000;      
gamma = 3;
gamma_c = 3;
d = 12;
d_min = 7;

t = 1;
TEST = cell(t,t);
mu = 0.9; %linspace(0.3,0.999,t);
epsilon = 0.15; % linspace(0.05,0.2,t);

for i = 1:t

    % GENERAZIONE DELLA RETE
    flag = true;
    tic
    while flag
    [A,~,c,dd] = LFR(n,d,mu(i),gamma, gamma_c, d_min);
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
    [opinionHistory, it] = HK(A, W, opin, confidence); % DIFFUSIONE
    toc
    % 
    TEST{i,j} = opinionHistory;
    % fprintf('mu = %d , eps = %d \n',mu(i),epsilon(j));

    end

    fprintf('\n');
end

%% Test dei grafici
op = TEST{1,1};
dinamicaOpinione(op,0)
% plotOpinionBins(op,0.002,1)
% heatmap2D(op)
