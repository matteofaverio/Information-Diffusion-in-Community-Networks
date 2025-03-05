function A = rewiring(A, c, mu, maxit, toll)

flag = 1;
n = size(A,1);
deg = sum(A,2);

% Calcola current_int per ogni nodo e l'errore
current_int = zeros(n,1);
for i = 1:n
    for j = [1:i-1,i+1:n]
        if (c(i) == c(j) && A(i,j) == 1)
            current_int(i) = current_int(i) + 1;
        end
    end
end

errors = current_int./deg - mu ;
flag
%% PER I NODI CON current_int > desired_int (err>0)
% Per ottimizzare la ricerca delle coppie di nodi che massimizzano 
%  l'errore, si implementano le seguenti strutture dati.

% Crea una cella di tante matrici quante sono le comunità. 
% Ogni matrice ha 2 colonne ( idx , err ) e tante righe quanti gli elementi
% della comunità.
comm_index = unique(c);
num_comm = length(comm_index);
err_comm = cell(1,num_comm);
for i = 1:n
    err_comm{c(i)} = [err_comm{c(i)};i,errors(i)];
end

% Crea una cella di tante matrici quante sono le comunità.
% Ogni matrice è la matrice di adiacenza della comunità
A_comm = cell(1,num_comm);
for i = 1:num_comm
    A_comm{i} = A(err_comm{i}(:,1),err_comm{i}(:,1));
end
%%
% Crea una matrice per le coppie i cui errori sommati sono massimi 
% comunità per comunità
% La matrice ha 3 righe: 
% ( err_max ; idx_max1 ; idx_max2 ) 
% e tante colonne quante sono le comunità.
ERR = zeros(3,num_comm);
for i = 1:num_comm
    ERR(1:3,i) = find_couple_max(err_comm{i},A_comm{i});
end
flag
%% PER I NODI CON current_int < desired_int (err<0)

% creo struttura dati ordinata 

%%
it = 0;

while( norm(errors)>toll && it<maxit )
    
    % Scegli le due coppie ottimali per errori positivi
    [~, idx] = maxk(ERR(1, :), 2);
    % Trova gli indici a,b,c,d dei 4 nodi corrispondenti. (a,b) e (c,d) 
    % sarranno coppie di nodi appartenenti a due comunità distinte
    a = ERR(2,idx(1)); b = ERR(3,idx(1)); c = ERR(2,idx(2)); d = ERR(3,idx(2));
    % Stacca a-b, stacca c-d, collega a-c, collega b-d per la matrice A e
    % per le matrici A_comm coinvolte
    A(a,b) = 0; A(b,a) = 0; A(c,d) = 0; A(d,c) = 0; 
    A(a,c) = 1; A(c,a) = 1; A(b,d) = 1; A(d,b) = 1;
    p = find( err_comm{idx(1)}(:,1) == a, 1); r = find( err_comm{idx(1)}(:,1) == b, 1);
    s = find( err_comm{idx(2)}(:,1) == c, 1); t = find( err_comm{idx(1)}(:,1) == d, 1);
    A_comm{idx(1)}(p,r) = 0; A_comm{idx(1)}(r,p) = 0;
    A_comm{idx(1)}(s,t) = 0; A_comm{idx(1)}(t,s) = 0;
    A_comm{idx(1)}(a,s) = 1; A_comm{idx(1)}(s,a) = 1;
    A_comm{idx(1)}(r,t) = 1; A_comm{idx(1)}(t,r) = 1;
    % Togliere 1 a current_int per a,b,c,d
    current_int(a) = current_int(a) - 1; current_int(b) = current_int(b) - 1; 
    current_int(c) = current_int(c) - 1; current_int(d) = current_int(d) - 1;
    % Ricalcolo errore per a,b,c,d
    err_comm{idx(1)}(p,2) = err_comm{idx(1)}(p,2) - 1;
    err_comm{idx(1)}(r,2) = err_comm{idx(1)}(r,2) - 1;
    err_comm{idx(2)}(s,2) = err_comm{idx(2)}(s,2) - 1;
    err_comm{idx(2)}(t,2) = err_comm{idx(2)}(t,2) - 1;


    
    

    errors = current_int./deg - mu ;
    it = it + 1
    
end