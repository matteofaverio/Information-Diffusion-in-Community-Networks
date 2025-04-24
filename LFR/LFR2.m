function [A,A_intra,c,dd] = LFR2(n,d,mu,gamma, gamma_c, d_min)
%%
% network_LFR creates a Lancichinetti-Fortunato-Radicchi network with n
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

% Generazione dei gradi dei nodi in base a una power law distribution
[dd, d_max] = powerLaw_degree(n,gamma,d_min,d);
%fprintf('d_max = %3d\n',d_max)

% Generazione delle dimensioni delle comunità
[S,N] = powerLaw_communities(n,d_min,d_max,gamma_c);
%fprintf('Numero di comunità generate: %d\n', N)

%% Assegnazione di ogni nodo a una comunità


A_intra = zeros(n); A = zeros(n);

[~, sorted_idx] = sort(dd, 'descend');
c = zeros(n, 1);
inhabits = zeros(N, 1);
S = S';

for ii = 1:n
    i = sorted_idx(ii);
    % Seleziona le comunità che hanno ancora spazio e che soddisfano S >= ceil(mu * dd(i))
    available_comms = find(inhabits < S);
    big_enough_comms = find(S >= ceil(mu * dd(i)));
    viable_comms = intersect(available_comms,big_enough_comms);
    if isempty(viable_comms)
        fprintf('Impossibile inserire i nodi nelle comunità generate\n')
        return
    end
    % Tra quelle valide, scegli quella con il maggiore spazio residuo
    [~, best_idx] = max(S(viable_comms) - inhabits(viable_comms));
    chosen_comm = viable_comms(best_idx);
    c(i) = chosen_comm;
    inhabits(chosen_comm) = inhabits(chosen_comm) + 1;
end
%fprintf('Assegnazione completata con successo.\n');

%% Creazione dei collegamenti tra i nodi all'interno delle comunità


for comm = 1:N
    nodes_comm = find(c == comm);

    x = dd(nodes_comm)*mu;
    f = (x - floor(x));

    target_internal =  floor(x) + ( f > (1-mu)^2 );

    % Costruzione della tabella per la comunità
    T = zeros(length(nodes_comm),2);
    for i = 1:length(nodes_comm)
        T(i,:) = [nodes_comm(i) , target_internal(i)];
    end

    t = T(T(:,2) > 0, :);
    how_many = size(t,1);

    flag = true; tried = 1;

    % Costruzione di un set di possibili nodi da sostituire ad eventuali
    % nodi problematici
    set = nodes_comm(dd(nodes_comm) > d_min);
    [~, idx_sort] = sort(dd(set),'descend'); 
    set = set(idx_sort);
    max_tried = length(set);
 
    while flag && how_many > 1

        [~, idx] = sort(t(:,2), 'descend');
        u = t(idx(1), 1);
        v = t(idx(2), 1);

        if A_intra(u, v) > 0
            success = false;
            count = 2;
            while ~success && count < how_many
                count = count + 1;
                candidate = t(idx(count), 1);
                if (A_intra(u,candidate) == 0)
                    success = true;
                end
            end
            if success
                A_intra(u, candidate) = 1;
                A_intra(candidate, u) = 1;
                t(idx(1),2) = t(idx(1),2) - 1;
                t(idx(count),2) = t(idx(count),2) - 1;
            else 
               t(idx(1),:) = [];
               if ~isempty(set)
                   while ~isempty(t(:,1)==set(tried)) && tried < max_tried
                          tried = tried + 1;
                          success = true;
                   end
                   if success
                       t(end+1,:) = [set(tried),1];
                   else
                       flag = false;
                   end
               end
            end

        else 
            A_intra(u, v) = 1;
            A_intra(v, u) = 1;
            t(idx(1),2) = t(idx(1),2) - 1;
            t(idx(2),2) = t(idx(2),2) - 1;
        end

        t = t(t(:,2) > 0, :);
        how_many = size(t,1);

    end
end
%fprintf('Collegamenti intra-comunità generati.\n');

%% creazione dei collegamenti fra le comunità

A_inter = zeros(n);
external_degree = dd - sum(A_intra, 2);

T = zeros(n,2);
for i = 1:n
    T(i,:) = [i , external_degree(i)];
end

t = T(T(:,2) > 0, :);
how_many = size(t,1);

while how_many > 1

    [~, idx] = sort(t(:,2), 'descend');
    u = t(idx(1), 1);
    v = t(idx(2), 1);

    if c(u) == c(v) || A_inter(u, v) > 0
        success = false;
        count = 2;
        while ~success && count < how_many
            count = count + 1;
            candidate = t(idx(count), 1);
            if (A_inter(u,candidate) == 0)
                success = true;
            end
        end
        if success
            A_inter(u, candidate) = 1;
            A_inter(candidate, u) = 1;
            t(idx(1),2) = t(idx(1),2) - 1;
            t(idx(count),2) = t(idx(count),2) - 1;
        else
            t(idx(1),:) = [];
        end
    else 
        A_inter(u, v) = 1;
        A_inter(v, u) = 1;
        t(idx(1),2) = t(idx(1),2) - 1;
        t(idx(2),2) = t(idx(2),2) - 1;
    end

    t = t(t(:,2) > 0, :);
    how_many = size(t,1);

end

%fprintf('Collegamenti inter-comunità generati.\n');


%% creazione matrice risultante

A = A_intra + A_inter;

