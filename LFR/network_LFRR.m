function [A,c,L,dd] = network_LFRR(n,d,mu,gamma, gamma_c, d_min)
%network_LFR creates a Lancichinetti-Fortunato-Radicchi network with n
%nodes,average connectivity d,and mixing parameter mu. The nodes
%distribution follows a power law with exponent gamma and the community
%scale follows a power law with exponent gamma_c
%   [A] = network_LFR(n,m)
%   Input values:
%   n: number of desired nodes
%   gamma:
%   gamma_c:
%   d_min:
% 
%   Output Values:
%   A: adjacency matrix of the generated graph
%   c: vector indicating which community each node belongs to

A = zeros(n);

% Generazione dei gradi dei nodi in base a una power law distribution
[dd, d_max] = powerLaw_degree(n,gamma,d_min,d);
fprintf('d_max = %3d\n',d_max)

% Generazione delle dimensioni delle comunità
[S,N] = powerLaw_communities(n,d_min,d_max,gamma_c);
fprintf('Numero di comunità generate: %d\n', N)

%% Creazione dei collegamenti tra i nodi
L = dd;

while sum(L) > 1

    first_pick = randi([1, sum(L)]);
    % trovo il nodo a cui corrisponde la prima scelta
    k = 0;
    while sum(L(1:k)) < first_pick
        k = k + 1;
    end
    L(k) = L(k)-1;
    
    flag = k;
    while( flag==k ) % esco dal ciclo solo se h è diverso da k (no self loops)
        second_pick = randi([1, sum(L)]);
        % trovo il nodo a cui corrisponde la seconda scelta
        h=0;
        while sum(L(1:h)) < second_pick
            h = h + 1;
        end
        flag = h;
    end

    L(h) = L(h) - 1;
    % creo il collegmento tra i due nodi
    A(k,h) = 1;
    
end
fprintf('Nodi collegati\n')

%% Assegnazione di ogni nodo a una comunità

% vettore contente la comunità a cui ogni nodo è assegnato
c = zeros(n,1);
% vettore contenente i nodi che devono ancora venire assegnati a una
% comunità
homeless = ones(n,1);
% vettore contenente il numero di nodi presenti in ogni comunità
inhabits = zeros(N,1);

% assegno ogni nodo a una comunità che può contenerlo
% ripeto finchè a tutti i nodi non è stata assegnata una comunità:
while sum(homeless) > 0
    % scelgo un nodo a caso senza una comunità
    to_pick = randi([1, sum(homeless)]);
    still_to_pick = find(homeless==1);
    i = still_to_pick(to_pick);

    % scelgo a random una comunità a cui associare il mio nodo
    comm = randi([1,N]);
    % assegno il mio nodo alla comunità solo se il numero di elementi di
    % questa supera il grado interno del nodo, ossia il numero di
    % collegamenti che esso stringe internamente alla comunità
    if S(comm) > mu*dd(i)
        c(i) = comm;
        homeless(i) = 0;
        % se la comunità è già piena sostituisco un elemento a caso con il 
        % nodo estratto
        if inhabits(comm) == S(comm)
            % caccio il nodo meno problematico (col grado minore)
            nodes_in_community = find(c == comm);
            [~ ,idx_to_kick] = min(dd(nodes_in_community));
            c(nodes_in_community(idx_to_kick)) = 0;
            homeless(nodes_in_community(idx_to_kick)) = 1;
        else
            inhabits(comm) = inhabits(comm) + 1;
        end
    end
end


A = A | A';
fprintf('Nodi assegnati alle comunità\n')

%% Reshuffle dei collegamenti 
% l'obiettivo è ottenere la giusta proporzione fra collegamenti con
% l'esterno e con l'interno della comunità per ogni nodo

A = rewiring(A,c,mu,10*n);

end

