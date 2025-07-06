%% Prova su una sola rete;
clc, clear, close all
folderName = 'test_mu050_processed';
data = load("C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Data\data_mu01.mat");

if ~isfolder(folderName)
    % tenta di risolvere folderName come folder in MATLAB path
    p = which(folderName);
    if isempty(p)
        error('Cartella %s non trovata né come directory né sul path.', folderName);
    else
        folderName = fileparts(p);  % prende la directory in cui risiede un .m con quel nome
    end
end

%% estrazione della matrice di adiecenza e vettore comunità
A = data.M;
A = A{1};
C = data.C;
C = C{1};


%% estrazione delle opinioni finali
opinions = load('opinions_mu050_eps001.mat');
op = opinions.outputCell{2,1};
op = double(op)/1000;

%%
[op_mean,op_mode,op_var] = opin_in_comm(C,op);
op_stdev = sqrt(op_var); 
[H] = shannonEntropy(C,op);
[BC] = bimodalityCoefficient(C,op);
CI = concentrationIndex(C,op);

%% Versione artigianale
clc, clear, close all
% folderName = 'test_mu050_processed';
results0 = community_breakdown('test_mu095_processed');

% Versione GPT

results1 = community_breakdown_fast('test_mu095_processed');

% Versione GEMINI

results2 = community_breakdown_opt('test_mu095_processed');

% Versione GEMINI 2

results3 = community_breakdown_definitivo('test_mu095_processed');

%%
% a = results0.ShannonEntropy{3};
% sum(a - results3.ShannonEntropy(:,3))
% sum(results0.ConcentrationIndex{3} - results3.ConcentrationIndex(:,3))
% sum(results0.BimodalityCoefficient{3} - results3.BimodalityCoefficient(:,3))
for i = 1:25
    a(i) = max(results0.ShannonEntropy{i} - results2.ShannonEntropy{i});
end
max(a)

%%
clc; clear; close all;
mus = 50:5:95;
nmus = numel(mus);
tic
for i = 1:nmus
    fprintf('===== inizio mu %03d =========================================\n', mus(i));
    mu = sprintf('%03d', mus(i));
    folderName = ['test_mu' mu '_processed'];

    % esegui l’analisi vera (qui è un placeholder)
    results = community_breakdown_definitivo(folderName);

    % path dove salvare
    folder = 'C:\Users\giogu\OneDrive - Politecnico di Milano\Desktop\Poli\Terzo anno\Tesi\Information-Diffusion-in-Community-Networks\Esperimento\Plot functions\community analysis\risultati analisi BC modificata';
    if ~exist(folder,'dir')
        mkdir(folder);
    end

    % costruisco nome e salvo
    nomeFile = ['analisi_comunita_Gem_mu' mu '.mat'];
    fname    = fullfile(folder, nomeFile);
    save(fname, 'results');

    fprintf('===== fine mu %03d ===========================================\n', mus(i));
end
toc

%%
clc, clear,close all
folderName = 'risultati analisi Gemini';

[H, BC, CI] = community_plots(folderName);
%%
plot_H_BC(H,BC)
%%
clc, clear
folderName = 'risultati analisi';
community_plots_artigianale(folderName)

