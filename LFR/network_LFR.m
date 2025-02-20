function [A,AA,c,homeless,L,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min)
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

% Generazione delle dimensioni delle comunità
[S,N] = powerLaw_communities(n,d_min,d_max,gamma_c);
fprintf('Numero di comunità generate: %d\n', N)

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
        % se la comunità è già sostisco un elemento a caso con il nodo
        % estratto
        if inhabits(comm) == S(comm)
            % scelgo elemento da cacciare a caso
            to_kick = randi([1, S(comm)]);
            nodes_in_community = find(c == comm);
            c(nodes_in_community(to_kick)) = 0;
            homeless(nodes_in_community(to_kick)) = 1;
        else
            inhabits(comm) = inhabits(comm)+1;
        end
    end
    % fprintf('nodo %d assegnato alla comunità %d\n', i,comm)
    %f = f+1
end

%% Creazione dei collegamenti tra i nodi all'interno delle comunità

A = zeros(n);

for i = 1:N
    % la lista L contiene il numero di collegamenti che devono ancora essere
    % formati dal nodo all'interno della sua comunità
    L = round(mu*dd.*(c == i));
    
    while sum(L) > 1
    
        first_pick = randi([1, sum(L)]);
        % trovo il nodo a cui corrisponde la prima scelta
        k = 0;
        while sum(L(1:k)) < first_pick
            k = k + 1;
        end
        L(k) = L(k)-1;
  
        % trovo il nodo a cui corrisponde la seconda scelta
        h = 0;

        second_pick = randi([1, sum(L)]);
        while sum(L(1:h)) < second_pick
            h = h + 1;
        end

        % while h ~=0 || h~=k
        %     second_pick = randi([1, sum(L)]);
        %     while sum(L(1:h)) < second_pick
        %         h = h + 1;
        %     end
        % end
        L(h) = L(h) - 1;
        % creo il collegmento tra i due nodi
        A(k,h) = A(k,h) + 1;
        % rimuovo i collegamenti effettuati dalla lista
    end
end
AA = A;
%% creazione dei collegamenti fra le comunità
L = round((1-mu)*dd);

while sum(L) > 1

    first_pick = randi([1, sum(L)]);
    % trovo il nodo a cui corrisponde la prima scelta
    k = 0;
    while sum(L(1:k)) < first_pick
        k = k + 1;
    end
    L(k) = L(k)-1;
    % scelgo il secondo elemento in modo che non faccia parte della stessa
    % comunità del primo
    L_out = L.*( 1-( c==c(k)));
    second_pick = randi([1, sum( L_out) ]);
    % trovo il nodo a cui corrisponde la seconda scelta
    h = 0;
    while sum(L_out(1:h)) < second_pick
        h = h + 1;
    end
    L(h) = L(h) - 1;
    % creo il collegmento tra i due nodi
    A(k,h) = A(k,h) + 1;
end
end
