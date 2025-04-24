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

%%
