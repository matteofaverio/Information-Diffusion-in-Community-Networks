function [H] = shannonEntropy_opt(c, op, nbins, normalise)
% Versione DEFINITIVA di shannonEntropy. Corregge l'errore di dimensione.
if nargin < 3 || isempty(nbins),      nbins     = 25;   end
if nargin < 4 || isempty(normalise),  normalise = true; end

% 1. VALIDAZIONE E FILTRAGGIO (già corretto in precedenza)
c = c(:);
op = op(:);
if numel(c) ~= numel(op)
    error('shannonEntropy_opt:InputSizeMismatch', ...
          'I vettori di input c (comunità) e op (opinioni) devono avere lo stesso numero di elementi.');
end
valid_idx = ~isnan(op) & (op >= 0) & (op <= 1);
op_valid = op(valid_idx);
c_valid  = c(valid_idx);

if isempty(c_valid)
    H = [];
    return;
end

% --- INIZIO FIX per errore horzcat ---
% 2. MAPPATURA GRUPPI: Usiamo il TERZO output di 'unique' per ottenere l'indice per ogni nodo.
%    'groups' conterrà i nomi unici dei gruppi (es. [1, 2, 5, 10])
%    'group_idx' conterrà un vettore della stessa dimensione di c_valid, con l'indice del gruppo (es. [1,1,2,3,3,4...])
[groups, ~, group_idx] = unique(c_valid, 'stable');
num_groups = numel(groups);
% --- FINE FIX ---

% 3. BINNING DELLE OPINIONI
edges   = linspace(0, 1, nbins + 1);
bin_idx = discretize(op_valid, edges);

% 4. ACCUMARRAY: Ora 'group_idx' e 'bin_idx' hanno la stessa dimensione e la concatenazione funziona.
counts = accumarray([group_idx, bin_idx], 1, [num_groups, nbins]);

% 5. CALCOLO ENTROPIA (invariato)
sum_counts = sum(counts, 2);
sum_counts(sum_counts == 0) = 1; 
p = counts ./ sum_counts;

log_p = log(p);
log_p(p == 0) = 0;
H = -sum(p .* log_p, 2);

if normalise
    H = H / log(nbins);
end
end