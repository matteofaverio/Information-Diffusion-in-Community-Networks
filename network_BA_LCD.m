function [A] = network_BA_LCD(n,m)
%network_BA_LCD creates a Barabàsi_Albert network with n nodes using
%linearized chord diagram
%   [A] = network_BA(n)
%   Input values:
%   n: number of desired nodes
%   m: half of the average degree
% 
%   Output Values:
%   A: adjacency matrix of the generated graph

AA = zeros(n*m);
AA(1,1) = 1;

% Costruzione di del grafo con m = 1
for i = 1:m*n

    newlink = randi([1,2*i-1]);

    % calcola il grado di ogni nodo
    d = 2*sum(AA,2);
    d(i) = 1;
    j = 1;
    while sum(d(1:j)) < newlink
            j = j+1;
    end

    % poichè ci sono elementi nulli mi riporto all'ultimo elemento
    % nullo che soddisfa la condizione di cui sopra
    while d(j) == 0
        j = j-1;
    end
    AA(i,j) = 1;

end
% copio la parte inferiore della matrice per renderla simmetrica
AA = AA+tril(AA)'-diag(diag(AA));

% Costruisco il grafico con m > 1 a partire da quello con m = 1
A = zeros(n);

for i = 1:n
    for j = 1:n
        A(i,j) = sum(sum(AA( (i - 1)*m + 1:i*m,(j - 1)*m + 1:j*m) ));
    end
end


end

