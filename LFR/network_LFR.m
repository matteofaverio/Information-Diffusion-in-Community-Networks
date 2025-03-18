function [A,AA,c,dd,s] = network_LFR(n,d,mu,gamma, gamma_c, d_min)
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
%fprintf('d_max = %3d\n',d_max)

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
    % fprintf('Nodo: %d, Grado: %d, Dimensione: %d\n',i, ceil(dd(to_pick)),S(comm))
    % assegno il mio nodo alla comunità solo se il numero di elementi di
    % questa supera il grado interno del nodo, ossia il numero di
    % collegamenti che esso stringe internamente alla comunità
    if S(comm) >= ceil(mu*dd(i))
        homeless(i) = 0;
        % se la comunità è già piena sostuisco un elemento a caso con il 
        % nodo estratto
        if inhabits(comm) == S(comm)
            % scelgo elemento da cacciare a caso
            to_kick = randi([1, S(comm)]);
            nodes_in_community = find(c == comm);
            c(nodes_in_community(to_kick)) = 0;
            homeless(nodes_in_community(to_kick)) = 1;
        else
            inhabits(comm) = inhabits(comm)+1;
        end
        c(i) = comm;
    end
    % fprintf('nodo %d assegnato alla comunità %d\n', i,comm)
end
%fprintf('Nodi assegnati alle comunità\n');

%% Creazione dei collegamenti tra i nodi all'interno delle comunità

% A = zeros(n);
% 
% for i = 1:N
%     % la lista L contiene il numero di collegamenti che devono ancora essere
%     % formati dal nodo all'interno della sua comunità
%     L = round(mu*dd.*(c == i));
% 
%     while sum(L) > 1
% 
%         first_pick = randi([1, sum(L)]);
%         % trovo il nodo a cui corrisponde la prima scelta
%         k = 0;
%         while sum(L(1:k)) < first_pick
%             k = k + 1;
%         end
%         L(k) = L(k)-1;
% 
%         % trovo il nodo a cui corrisponde la seconda scelta
%         h = 0;
% 
%         second_pick = randi([1, sum(L)]);
%         while sum(L(1:h)) < second_pick
%             h = h + 1;
%         end
% 
%         % while h ~=0 || h~=k
%         %     second_pick = randi([1, sum(L)]);
%         %     while sum(L(1:h)) < second_pick
%         %         h = h + 1;
%         %     end
%         % end
%         L(h) = L(h) - 1;
%         % creo il collegmento tra i due nodi
%         A(k,h) = A(k,h) + 1;
%         % rimuovo i collegamenti effettuati dalla lista
%     end
%     disp(sum(L))
% end
%%
residual_links = zeros(n,1);
s = 0;
L1 = dd;
for i = 1:N

    deg = round( (mu)*dd.*(c == i));
    L1 = L1-deg;
    s = s + sum(deg);
    degrees = unique(sort(deg,'descend'));

    L = zeros(n,1); 

    for j =  1:( length(degrees)-1 )
        d = degrees(j);
        L = L + deg.*(deg == d);

        if sum(L > 0) < 2
            continue;
        end

        while min(L(L>0)) > degrees(j+1)

                first_pick = randi([1, sum(L)]);
                % trovo il nodo a cui corrisponde la prima scelta
                k = find(cumsum(L) >= first_pick, 1);
                L(k) = L(k)-1;

                % faccio in modo di scegliere un nodo diverso da quello
                % precedente
                L_mod = L;
                L_mod(k) = 0; 

                second_pick = randi([1, sum(L_mod)]);
                % trovo il nodo a cui corrisponde la seconda scelta
                h = find(cumsum(L_mod) >= second_pick, 1);

                % rimuovo i collegamenti effettuati dalla lista
                L(h) = L(h) - 1;

                % creo il collegmento tra i due nodi
                A(k,h) = A(k,h) + 1;
        end
    end
        L = L + deg.*(deg == degrees(end));
    while sum(L) > 1
            first_pick = randi([1, sum(L)]);
            % trovo il nodo a cui corrisponde la prima scelta
            k = find(cumsum(L) >= first_pick, 1);
            L(k) = L(k)-1;

            % faccio in modo di scegliere un nodo diverso da quello
            % precedente
            L_mod = L;
            L_mod(k) = 0; 

            second_pick = randi([1, sum(L_mod)]);
            % trovo il nodo a cui corrisponde la seconda scelta
            h = find(cumsum(L_mod) >= second_pick, 1);

            % rimuovo i collegamenti effettuati dalla lista
            L(h) = L(h) - 1;

            % creo il collegmento tra i due nodi
            A(k,h) = A(k,h) + 1;
    end
            residual_links = residual_links + L;
end
AA = A;
fprintf('Comunità create\n')
%% creazione dei collegamenti fra le comunità
% L = round((1-mu)*dd) + residual_links;
L1 = L1 + residual_links;
s = s+sum(L1);

while sum(L) > 1

    first_pick = randi([1, sum(L)]);
    % trovo il nodo a cui corrisponde la prima scelta
    k = find(cumsum(L) >= first_pick, 1);
    % scelgo il secondo elemento in modo che non faccia parte della stessa
    % comunità del primo
    L_out = L;
    L_out(c == c(k)) = 0;
    if sum(L_out) > 0
        second_pick = randi([1, sum( L_out) ]);
    
        % trovo il nodo a cui corrisponde la seconda scelta
        h = find(cumsum(L_out) >= second_pick);
        L(h) = L(h) - 1;
        L(k) = L(k) - 1;
        % creo il collegmento tra i due nodi
        A(k,h) = A(k,h) + 1;
    else
        L(k) = L(k) + 1;
    end
end
%fprintf('Comunità collegate\n')

end
