function  results = community_breakdown(folderName)
%   Serve a elaborare la distribuzione dell'opinione nelle comunità per il
%   set di esperimenti conenuti in foldername. carica tutti i .mat in folderName,
%   estrae dalle 200 matrici in ciascun file l'ultima colonna (stato finale)
%   sulla base del risultato calcola media, valore modale (su 25 bin), e
%   varianza normalizzata
tic
% caricamento dati
epsilons   = 0.01:0.01:0.25;      % 25 valori di epsilon
nEps       = numel(epsilons);

% inizializzazione variabili output
entropy_data = cell(nEps);
CI_data = cell(nEps);
BC_data = cell(nEps);

fprintf('--- Caricamento files ---\n')
% recupero dati sulle comunità
direction = 'C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Data';
% Estrai le tre cifre dopo "mu"
tok = regexp(folderName, 'test_mu(\d{3})_processed', 'tokens', 'once');
% Convertilo in numero intero
val = str2double(tok{1});
% Trasforma in indice ordinale
idx = (val - 50)/5 + 1;
idxStr = sprintf('%02d', idx);
% nome del file da caricare
fileName = fullfile(direction, "data_mu" + idxStr + ".mat");

data = load(fileName);

% estraggo i dati sulle comunità
Communities = data.C;

% --- Elenco e ordinamento file ---
files = dir(fullfile(folderName,'*.mat'));
%files = files(~strncmp({files.name}, '._', 2));
if numel(files)~=nEps
    error('Trovati %d file, ma ci aspettavamo %d (numero di epsilon)', ...
          numel(files), nEps);
end
names = sort({files.name});
fprintf('--- files caricati ---\n')
for i = 1:nEps
    entropy_data{i} = zeros(200,1); 
    BC_data{i} = zeros(200,1); 
    CI_data{i} = zeros(200,1);
    S = load(fullfile(folderName, names{i}));
    fn = fieldnames(S);
    ops  = S.(fn{1}); % cell array 3x200

    for j = 1:200
        % estrazione delle opinioni finali e riconversione in double
        col  = double(ops{2,j})/1000;  

        % calcolo dei parametri medi
        c = Communities{j};
        entropy_data{i}(j) = mean(shannonEntropy(c,col));
        BC_data{i}(j) = mean(bimodalityCoefficient(c,col));
        CI_data{i}(j) = mean(concentrationIndex(c,col));

    end
    fprintf('--- epsilon %.2f completato ---\n', epsilons(i))
end
results = struct('ShannonEntropy', {entropy_data}, ...
                 'ConcentrationIndex', {CI_data}, ...
                 'BimodalityCoefficient', {BC_data});
toc
