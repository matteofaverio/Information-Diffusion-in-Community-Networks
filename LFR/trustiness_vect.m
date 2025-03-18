function [W] = trustiness_vect(A,c,foreign_trust)
%TRUSTINESS cretaes the trustiness matrix for an undirected graph.
% the trust of j in i is assigned randomly witch probability of i to
% influence j proprtional to the ratio between the degrees of i and j

n = size(A,1);           % Numero di righe (e colonne) di A
W = zeros(n);            % Inizializza la matrice W a zeri
deg = sum(A, 2);         % Calcola il vettore dei gradi: somma degli elementi per ogni riga

% Ottieni gli indici (i,j) per le coppie dove i > j (parte triangolare inferiore)
[I, J] = find(tril(ones(n), -1));

% Numero totale di coppie (i,j) con i > j
N = numel(I);

% Genera un vettore di N numeri casuali (uno per ciascuna coppia)
randVals = rand(N, 1);

% Per ogni coppia, calcola il "total" che è la somma dei gradi dei due nodi:
% total = deg(i) + deg(j)
total = deg(I) + deg(J);

% Calcola w per ciascuna coppia:
% w = (total * random_value - deg(i)) / total
% Qui, ".*" e "./" indicano moltiplicazione e divisione elemento per elemento.
w = (total .* randVals - deg(I));

% Ora, per ogni coppia abbiamo un valore w:
% - Se w > 0, allora secondo la regola vogliamo assegnare W(i,j) = w.
% - Se w <= 0, allora vogliamo assegnare W(j,i) = -w.

% Trova gli indici delle coppie in cui w è positivo (w > 0)
idx_pos = w > 0;
% Gli altri saranno quelli in cui w è negativo o zero
idx_neg = ~idx_pos;

% Utilizziamo "sub2ind" per convertire gli indici (riga, colonna) in indici lineari
% e assegnare i valori a W.
% Per le coppie con w positivo:
W(sub2ind([n, n], I(idx_pos), J(idx_pos))) = w(idx_pos);
% Per le coppie con w negativo o zero:
W(sub2ind([n, n], J(idx_neg), I(idx_neg))) = -w(idx_neg);

W = W.*A;


% maschera booleana, M(i,j) == 1 se i e j fanno parte della stessa comunità
M = (c == c');
% maschera coi nodi che non fanno parte della stessa comunità
M1 = 1-M;
M1 = foreign_trust*M1;
M = M+M1;

W = W.*M;

% Normalizzo la matrice di trustiness in modo che la somma dei lati
% entranti in ogni nodo sia 1
colSum = sum(W,1);
W = W./colSum;
W(isnan(W)) = 0;

end

