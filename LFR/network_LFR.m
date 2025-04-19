function [A,AA,c,dd] = network_LFR(n,d,mu,gamma, gamma_c, d_min)
%network_LFR creates a Lancichinetti-Fortunato-Radicchi network with n
%nodes,average connectivity d,and mixing parameter mu. The nodes
%distribution follows a power law with exponent gamma and the community
%scale follows a power law with exponent gamma_c
%   [A,AA,c,dd] = network_LFR(n,m)
%   Input values:
%   n: number of desired nodes
%   gamma:
%   gamma_c:
%   d_min:
% 
%   Output Values:
%   A: adjacency matrix of the generated graph
%   c: vector indicating which community each node belongs to

maxit = 1000;

A = zeros(n);
% Generazione dei gradi dei nodi in base a una power law distribution
[dd, d_max] = powerLaw_degree(n,gamma,d_min,d);
% fprintf('d_max = %3d\n',d_max)

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
it = 0;
while sum(homeless) > 0 && it < 10*n
    it = it+1;

    % scelgo un nodo a caso senza una comunità
    to_pick = randi([1, sum(homeless)]);
    still_to_pick = find(homeless==1);
    i = still_to_pick(to_pick);

    % scelgo a random una comunità a cui associare il mio nodo
    comm = randi([1,N]);

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
    
end
fprintf('numero tentativi assegnazione %d\n', it)
if it == n*10
    fprintf('Comunità non create')
    A = nan(n);
    AA = A;
    return
end
%fprintf('Nodi assegnati alle comunità\n');

%% Creazione dei collegamenti tra i nodi all'interno delle comunità

A = zeros(n);

% vettore che contiene tutti i collegamenti scartati
residual_links = zeros(n,1);
% vesttore contenente il numero di collegamenti intercomunitari
L1 = dd;

for i = 1:N
    % vettore contenente i collegamenti intracomunitari di ogni nodo
    deg = round( mu*dd.*(c == i));

    % rimuovo i collegamenti che verranno effettuati dai collegamenti
    % intracomunitari
    L1 = L1-deg;

    % lista ordinata dei numeri di collegamenti fra i nodi della comunità
    degrees = unique(sort(deg,'descend'));

    % lista che conterrà il numero di collegamenti mancanti per ogni nodo
    L = zeros(n,1); 
    
    % faccio un ciclo sui valori dei gradi interni alla comunità
    for j =  1:( length(degrees)-1 )
        d = degrees(j);

        % aggiungo i nodi con il grado corrente alla lista contenente i
        % collegamenti mancanti
        L = L + deg.*(deg == d);
        
        it = 0;     % contatore delle iterazioni

        % se la lista conitene un unico nodo passo al grado successivo
        if sum(L > 0) < 2
            continue;
        end
        

        while min(L(L>0)) > degrees(j+1) && it < maxit
                % scelgo il primo nodo da collegare
                first_pick = randi([1, sum(L)]);

                % trovo il nodo a cui corrisponde la prima scelta
                k = find(cumsum(L) >= first_pick, 1);

                % faccio in modo di scegliere un nodo diverso da quello
                % precedente
                L_mod = L;
                L_mod(k) = 0; 
                
                % se non ci sono altri nodi da collegare ripeto la prima
                % scelta
                if sum(L_mod) == 0
                    it = it+1;
                    continue
                end

                second_pick = randi([1, sum(L_mod)]);

                % trovo il nodo a cui corrisponde la seconda scelta
                h = find(cumsum(L_mod) >= second_pick, 1);
                   
                % se il collegamento esiste già ripeto il processo
                if A(k,h) + A(h,k) > 0
                    it = it+1;
                    continue
                end

                % rimuovo i collegamenti effettuati dalla lista
                L(h) = L(h) - 1;
                L(k) = L(k) - 1;

                % creo il collegmento tra i due nodi
                A(k,h) = A(k,h) + 1;
        end
    end

    it = 0;

    if it == maxit
        fprintf('Comunità non create')
        A = nan(n);
        return
    end

    % ripeto il procedimento per il grado più basso
    L = L + deg.*(deg == degrees(end));
    
    while sum(L) > 1 && it < maxit

            first_pick = randi([1, sum(L)]);
            % trovo il nodo a cui corrisponde la prima scelta
            k = find(cumsum(L) >= first_pick, 1);
          
             % faccio in modo di scegliere un nodo diverso da quello
            % precedente
            L_mod = L;
            L_mod(k) = 0; 
            
            % se non ci sono altri nodi da collegare ripeto la prima
            % scelta
            if sum(L_mod) == 0
                it = it+1;
                continue
            end

            second_pick = randi([1, sum(L_mod)]);

            % trovo il nodo a cui corrisponde la seconda scelta
            h = find(cumsum(L_mod) >= second_pick, 1);
            
            % se il collegamento esiste già ripeto il processo
            if A(k,h) + A(h,k) > 0
                it = it+1;
                continue
            end

            % rimuovo i collegamenti effettuati dalla lista
            L(h) = L(h) - 1;
            L(k) = L(k) - 1;

            % creo il collegmento tra i due nodi
            A(k,h) = A(k,h) + 1;
    end
    % resgistro in una lista tutti i collegamenti che non sono stati
    % effettuati
    residual_links = residual_links + L;
end
AA = A;
if it == maxit
    fprintf('Comunità non create')
    A = nan(n);
    AA = A;
    return
else
    fprintf('Comunità create\n')
end
%% creazione dei collegamenti fra le comunità
L1 = L1 + residual_links;

% contatore delle iterazioni
it = 0;

while sum(L1) > 1 && it < maxit

    first_pick = randi([1, sum(L1)]);

    % trovo il nodo a cui corrisponde la prima scelta
    k = find(cumsum(L1) >= first_pick, 1);

    % scelgo il secondo elemento in modo che non faccia parte della stessa
    % comunità del primo
    L_out = L1;
    L_out(c == c(k)) = 0;

    if sum(L_out) == 0
            it = it+1;
            continue
    end
    second_pick = randi([1, sum( L_out) ]);

    % trovo il nodo a cui corrisponde la seconda scelta
    h = find(cumsum(L_out) >= second_pick,1);

    if A(k,h) + A(h,k) > 0
            it = it+1;
            continue
    end

    L1(h) = L1(h) - 1;
    L1(k) = L1(k) - 1;

    % creo il collegmento tra i due nodi
    A(k,h) = A(k,h) + 1;
end

if it == maxit
    fprintf('Comunità non collegate')
    A = nan(n);
    return
end
%fprintf('Comunità collegate\n')
A = A+A';
end
