function [A_rewired,mu_new] = rewiringLFR(A, comm, mu, max_iter)

% A         : matrice di adiacenza simmetrica (non orientata)
% comm      : vettore di assegnazione comunità per ogni nodo (lunghezza N)
% mu        : frazione desiderata di collegamenti interni per ogni nodo
% max_iter  : numero massimo di iterazioni di rewiring

A = A - diag(diag(A));
N = size(A, 1);           
d = sum(A, 2);            

k_in_desired = mu * d;
k_in_current = zeros(N, 1);
for i = 1:N
    neighbors = find(A(i, :));
    k_in_current(i) = sum(comm(neighbors) == comm(i));
end

err = (k_in_desired - k_in_current);

[I, J] = find(triu(A));
internal_edges = [];
for idx = 1:length(I)
    i = I(idx); j = J(idx);
    if comm(i) == comm(j)
        internal_edges = [internal_edges; i, j, err(i)+err(j)]; 
    end
end

M = (comm == comm');
sameCommCounts = sum(A .* M, 2);
fractions = sameCommCounts ./ d;
mu_new = mean(fractions);

mu_error = mu_new-mu;

it_min = mu_error*20000;

iter = 0;
improved2 = true;

while iter < max_iter && improved2 && mu_error>1e-4

    num_int = size(internal_edges, 1);
    
    to_check_min = ones(num_int,1);

    improved2 = false;
    
    while ~improved2
        
        to_adjust = [];
        flag = 0;

        while ~flag
            [val_err_min,idx1] = min(internal_edges(:,3));
            flag = to_check_min(idx1);
            if val_err_min >= 0
                break;
            end
            if ~flag 
                internal_edges(idx1,3) = 0;
                to_adjust = [to_adjust;idx1];
            end
        end
        for i = to_adjust
            internal_edges(i,3) = err(internal_edges(i,1)) + err(internal_edges(i,2));
        end
        if val_err_min >= 0
            break;
        end 
        
        i = internal_edges(idx1,1); j = internal_edges(idx1,2);
        
        found2 = false;
        tried = 0;  
        to_check = ones(num_int,1);
        to_check(idx1) = 0;

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
            elseif internal_edges(idx1,3) + internal_edges(idx1,3) >= -1
                to_check(idx2) = 0;
            else
                if (A(i,h) || A(j,k)) found2 = true; idx = 2; 
                elseif (A(i,k) || A(j,h)) found2 = true; idx = 1;
                else found2 = true; idx = randi([1,2]);          
                end
            end
    
            if (tried == num_int-1)
                break; 
            end
    
        end
        
        if found2 
            improved2 = true;
            
            A(i,j) = 0; A(j,i) = 0; A(h,k) = 0; A(k,h) = 0;
            k_in_current([i,j,h,k]) = k_in_current([i,j,h,k]) - 1;
            err([i,j,h,k]) = err([i,j,h,k]) + 1;
            internal_edges([idx1,idx2],:) = [];
    
            if(idx==1)
                A(i,h) = 1; A(h,i) = 1; A(j,k) = 1; A(k,j) = 1;
    
            elseif(idx==2)
                A(i,k) = 1; A(k,i) = 1; A(j,h) = 1; A(h,j) = 1;
            end
        end
        
        if ~improved2 
            to_check_min(idx1) = 0;
        end
        
    end

    iter = iter + 1;

    if iter > it_min
        sameCommCounts = sum(A .* M, 2);
        fractions = sameCommCounts ./ d;
        mu_new = mean(fractions);
    
        mu_error = mu_new-mu;
    end

end

A_rewired = A;



end


