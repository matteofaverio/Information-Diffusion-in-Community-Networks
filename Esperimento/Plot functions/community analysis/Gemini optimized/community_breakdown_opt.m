function results = community_breakdown_opt(folderName)
% Versione ottimizzata di community_breakdown.

% --- INIZIO OTTIMIZZAZIONE ---
% Commento: Ho rimosso il percorso hardcoded. È meglio passarlo come argomento 
% o usare percorsi relativi. Qui lo definisco per compatibilità.
direction = 'C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Data';
% --- FINE OTTIMIZZAZIONE ---

tic
epsilons   = 0.01:0.01:0.25;
nEps       = numel(epsilons);

% Inizializzazione variabili output
entropy_data = cell(1, nEps);
CI_data      = cell(1, nEps);
BC_data      = cell(1, nEps);

fprintf('--- Caricamento dati e setup ---\n')
% Estrai l'indice mu dal nome della cartella
tok = regexp(folderName, 'test_mu(\d{3})_processed', 'tokens', 'once');
if isempty(tok)
    error('Il nome della cartella non corrisponde al pattern atteso: test_muXXX_processed');
end
val = str2double(tok{1});
idx = (val - 50)/5 + 1;
idxStr = sprintf('%02d', idx);

% Caricamento dati comunità (una sola volta)
fileName = fullfile(direction, ['data_mu' idxStr '.mat']);
data = load(fileName);
Communities = data.C;

% Elenco e ordinamento file
files = dir(fullfile(folderName,'*.mat'));
if numel(files) ~= nEps
    error('Trovati %d file, ma ci aspettavamo %d (numero di epsilon)', ...
          numel(files), nEps);
end
names = sort({files.name});
fprintf('--- Setup completato ---\n')

% --- INIZIO OTTIMIZZAZIONE ---
% CICLO PRINCIPALE: Se possiedi il Parallel Computing Toolbox, puoi
% sostituire 'for' con 'parfor' per un'enorme accelerazione, dato che ogni
% iterazione (ogni epsilon) è indipendente.
% parfor i = 1:nEps
for i = 1:nEps
    % Pre-allocazione per la singola iterazione di parfor
    entropy_i = zeros(200,1); 
    BC_i      = zeros(200,1); 
    CI_i      = zeros(200,1);

    S = load(fullfile(folderName, names{i}));
    fn = fieldnames(S);
    ops  = S.(fn{1}); % cell array 3x200
    
    for j = 1:200
        col = double(ops{2,j})/1000;
        c   = Communities{j};
        
        % Chiamata alle nuove funzioni vettorializzate
        entropy_i(j) = mean(shannonEntropy_opt(c, col));
        BC_i(j)      = mean(bimodalityCoefficient_opt(c, col));
        CI_i(j)      = mean(concentrationIndex_opt(c, col));
    end
    
    % Assegnazione dei risultati alla cella
    entropy_data{i} = entropy_i;
    BC_data{i}      = BC_i;
    CI_data{i}      = CI_i;

    fprintf('--- epsilon %.2f completato ---\n', epsilons(i))
end
% --- FINE OTTIMIZZAZIONE ---

results = struct('ShannonEntropy', {entropy_data}, ...
                 'ConcentrationIndex', {CI_data}, ...
                 'BimodalityCoefficient', {BC_data});
toc
end