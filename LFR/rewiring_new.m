function [A_rewired,err_agg] = rewiring_new(A, comm, mu, max_iter, tol)

% A         : matrice di adiacenza simmetrica (non orientata)
% comm      : vettore di assegnazione comunità per ogni nodo (lunghezza N)
% mu        : frazione desiderata di collegamenti interni per ogni nodo
% max_iter  : numero massimo di iterazioni di rewiring
% tol       : tolleranza sul deficit di collegamenti interni

N = size(A, 1);           % numero di nodi
d = sum(A, 2);            % vettore dei gradi

% Calcola il numero desiderato di collegamenti interni per ciascun nodo
k_in_desired = mu * d;

% Calcola il numero corrente di collegamenti interni per ciascun nodo
k_in_current = zeros(N, 1);
for i = 1:N
    neighbors = find(A(i, :));
    k_in_current(i) = sum(comm(neighbors) == comm(i));
end

% Calcola la discrepanza per ogni nodo
err = (k_in_desired - k_in_current);
err_agg = err;

% Estrai la lista degli archi esterni
% Si considerano solo i lati con i < j per evitare duplicazioni
[I, J] = find(triu(A));
external_edges = [];
for idx = 1:length(I)
    i = I(idx); j = J(idx);
    if comm(i) ~= comm(j)
        external_edges = [external_edges; i, j, err(i)+err(j)];
    end
end

% Estrai la lista degli archi interni
% Si considerano solo i lati con i < j per evitare duplicazioni
internal_edges = [];
for idx = 1:length(I)
    i = I(idx); j = J(idx);
    if comm(i) == comm(j)
        internal_edges = [internal_edges; i, j, err(i)+err(j)]; 
    end
end

iter = 0;
improved = true;
found2 = 1;

while iter < max_iter && improved

    improved = false;

    num_ext = size(external_edges, 1);
    
    to_check = ones(num_ext,1);
    [~,idx1] = max(external_edges(:,3));
    i = external_edges(idx1,1); j = external_edges(idx1,2);
    to_check(idx1) = 0;
    
    found1 = false;
    tried = 0;  

    while ~found1

        flag = 0;
        while ~flag
            idx2 = randi([1,num_ext]);
            flag = to_check(idx2);
        end
        h = external_edges(idx2,1);
        k = external_edges(idx2,2);
        tried = tried + 1;

        % controllo
        if (i==h || j==k || i==k || j == h) 
            to_check(idx2) = 0;
        elseif A(i,h) || A(j,k) || A(i,k) || A(j,h)
            to_check(idx2) = 0;
        elseif length(unique([comm(i),comm(j),comm(h),comm(k)]))==4
            to_check(idx2) = 0;
        else
            % calcolo l'effetto del caso1
            cc1 = [ (comm(i) == comm(h)) ; (comm(j) == comm(k)) ];
            d1 = k_in_current([i,j,h,k]) + [cc1;cc1];
            e1 = k_in_desired([i,j,h,k]) - d1;
            e_1 = sum(e1);
            % calcolo l'effetto del caso2
            cc2 = [ (comm(i) == comm(k)) ; (comm(j) == comm(h)) ];
            d2 = k_in_current([i,j,k,h]) + [cc2;cc2];
            e2 = k_in_desired([i,j,k,h]) - d2;
            e_2 = sum(e2);
            % cerco il caso con l'effetto massimo (errore minimo)
            [e_fin,idx] = min([e_1,e_2]);
            % se l'errore è ancora positivo, found = true
            if e_fin > 0
                found1 = true;
            else 
                to_check(idx2) = 0;
            end
        end
        if(tried == num_ext-1)
            break;
        end
    end
    
    if found1
        A(i,j) = 0; A(j,i) = 0; A(h,k) = 0; A(k,h) = 0;
        external_edges([idx1,idx2],:) = [];
        if(idx==1)
            
            A(i,h) = 1; A(h,i) = 1; A(j,k) = 1; A(k,j) = 1;
            k_in_current([i,j,h,k]) = d1;
            err([i,j,h,k]) = e1;
            if cc1(1) internal_edges(end+1,:) = [i,h,err(i)+err(h)]; end
            if cc1(2) internal_edges(end+1,:) = [j,k,err(j)+err(k)]; end
            if ~cc1(1) external_edges(end+1,:) = [i,h,err(i)+err(h)]; end
            if ~cc1(2) external_edges(end+1,:) = [j,k,err(j)+err(k)]; end

        elseif(idx==2)
            A(i,k) = 1; A(k,i) = 1; A(j,h) = 1; A(h,j) = 1;
            k_in_current([i,j,k,h]) = d2;
            err([i,j,k,h]) = e2;
            if cc2(1) internal_edges(end+1,:) = [i,k,err(i)+err(k)]; end
            if cc2(2) internal_edges(end+1,:) = [j,h,err(j)+err(h)]; end
            if ~cc2(1) external_edges(end+1,:) = [i,k,err(i)+err(k)]; end
            if ~cc2(2) external_edges(end+1,:) = [j,h,err(j)+err(h)]; end

        else
            fprintf('errore');
        end
    
    end
    
    num_int = size(internal_edges, 1);
    
    to_check = ones(num_int,1);
    [~,idx1] = min(internal_edges(:,3));
    i = internal_edges(idx1,1); j = internal_edges(idx1,2);
    to_check(idx1) = 0;
    
    found2 = false;
    tried = 0;  

    while ~found2

        flag = 0;
        while ~flag
            idx2 = randi([1,num_int]);
            flag = to_check(idx2);
        end
        h = internal_edges(idx2,1);
        k = internal_edges(idx2,2);
        tried = tried + 1;

        if comm(i) == comm(h)
            to_check(idx2) = 0;
        elseif (A(i,h) || A(j,k)) && (A(i,k) || A(j,h))
            to_check(idx2) = 0;
        elseif internal_edges(idx1,3) + internal_edges(idx1,3) > -4
            to_check(idx2) = 0;
        else
            if (A(i,h) || A(j,k)) found2 = true; idx = 2; 
            elseif (A(i,k) || A(j,h)) found2 = true; idx = 1;
            else found2 = true; idx = randi([1,2]);          
            end
        end
        if(tried == num_int-1)
            break;
        end
    end
    
    if found2 

        A(i,j) = 0; A(j,i) = 0; A(h,k) = 0; A(k,h) = 0;
        k_in_current([i,j,h,k]) = k_in_current([i,j,h,k]) - 1;
        err([i,j,h,k]) = err([i,j,h,k]) + 1;
        internal_edges([idx1,idx2],:) = [];

        if(idx==1)
            A(i,h) = 1; A(h,i) = 1; A(j,k) = 1; A(k,j) = 1;
            external_edges(end+1,:) = [i,h,err(i)+err(h)];
            external_edges(end+1,:) = [j,k,err(j)+err(k)];

        elseif(idx==2)
            A(i,k) = 1; A(k,i) = 1; A(j,h) = 1; A(h,j) = 1;
            external_edges(end+1,:) = [i,k,err(i)+err(k)];
            external_edges(end+1,:) = [j,h,err(j)+err(h)];

        else
            fprintf('errore');
        end
    
    end

    err_agg = [err_agg , err];
    improved = found1 || found2;
    iter = iter + 1;

end

A_rewired = A;
iter
end


