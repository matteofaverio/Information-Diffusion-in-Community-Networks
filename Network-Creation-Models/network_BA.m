function [A] = network_BA(n,m)
%network_BA creates a Barabàsi_Albert network with n with average degree
%2*m
%   [A] = network_BA(n,m)
%   Input values:
%   n: number of desired nodes
%   m: half of the average degree
% 
%   Output Values:
%   A: adjacency matrix of the generated graph

A = zeros(n);

% creo un nucleo iniziale di collegamenti
n0 = ceil(log2(n));
A(1:n0,1:n0) = ones(n0)-binornd(1,0.3*ones(n0),n0);

% inizializzo il numero di nodi
i = n0;

% aggiungo i nuovi agenti
while i < n

    i = i+1;

    % calcolo il grado di ogni agente
    k = sum(A,2);
    % calcolo il grado totale
    ktot = sum(k);
    % inizializzo il numero di collegamenti del nuovo nodo
    l = 0;
    
    % creo un vettore che servirà per estrarre il nuovo collegamento,
    % il vettore contiene il numero dei collegamento dei nodi a cui il
    % nuovo nodo non è collegato
    pool = k;
    pool(i) = 1;

    while l < m

        l = l+1;
        % genero un numero casuale che decide il nuovo collegamento, il
        % numero è generato da un'uniforme da 1 al numero totale di
        % collegamenti
        newlink = randi(sum(pool));
        
        j = 1;
        
        % trovo l'elemento a cui corrisponde il numero estratto 
        while sum(pool(1:j)) < newlink
            j = j+1;
        end

        % poichè ci sono elementi nulli mi riporto all'ultimo elemento
        % nullo che soddisfa la condizione di cui sopra
        while pool(j) == 0
            j = j-1;
        end
        A(i,j) = 1;
        
    end
end
A = tril(A) + tril(A)';
end


