function [A] = network_SM(n,k,p)
%network_SM creates a Small-World network with n with average degree
%2*k
%   [A] = network_SM(n,k,p)
%   Input values:
%   n: number of desired nodes
%   k: half of the average degree
%   p: probability that an edge is reshuffled
% 
%   Output Values:
%   A: adjacency matrix of the generated graph

% Genero la matrice iniziale in cui ogni nodo è connesso ai k nodi
% precedenti e ai k nodi successivi
A = zeros(n);

for i = 1:k
    A = A + diag(ones(1,n-i),i);
end

A(1:k,n-k+1:n) = A(1:k,n-k+1:n) + triu(ones(k));

A = A+triu(A);

% Opero il reshuffle dei lati 
B = zeros(n);

for i = 1:n
    for j = i:n 
        % parametro che determina se un lato connesso al lato i cambia
        % destinazione
        flag = binornd(1,p);

        if and(A(i,j) == 1 , flag == 1)
            B(i,j) = - 1;
            % scelta del nuovo vertice del lato
            r = randi(1,n);
            B(i,r) =  1;
        end
    end

end    

% applico la trasformazione alla matrice A
A = A+B;
A = triu(A)+triu(A)';
% normalizzo la matrice per avere solo elementi unitari
A = A./A;
A(isnan(A)) = 0;
% elimino gli elementi dalla diagonale principale
A = A-diag(diag(A));

end

