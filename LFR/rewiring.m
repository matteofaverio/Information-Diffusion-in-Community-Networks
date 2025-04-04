function [A_rewired,err_agg,var_agg] = rewiring(A, comm, mu, max_iter)

% A         : matrice di adiacenza simmetrica (non orientata)
% comm      : vettore di assegnazione comunità per ogni nodo (lunghezza N)
% mu        : frazione desiderata di collegamenti interni per ogni nodo
% max_iter  : numero massimo di iterazioni di rewiring

N = size(A, 1);           
d = sum(A, 2);            

k_in_desired = mu * d;
k_in_current = zeros(N, 1);
for i = 1:N
    neighbors = find(A(i, :));
    k_in_current(i) = sum(comm(neighbors) == comm(i));
end

err = (k_in_desired - k_in_current);
err_agg = mean(err);
var_agg = (sum((err).^2))./N;

% lista degli archi esterni
[I, J] = find(triu(A));
external_edges = [];
for idx = 1:length(I)
    i = I(idx); j = J(idx);
    if comm(i) ~= comm(j)
        external_edges = [external_edges; i, j, err(i)+err(j)];
    end
end

% lista degli archi interni
internal_edges = [];
for idx = 1:length(I)
    i = I(idx); j = J(idx);
    if comm(i) == comm(j)
        internal_edges = [internal_edges; i, j, err(i)+err(j)]; 
    end
end

iter = 0;
improved = true;
count1 = 0;
count2 = 0;

while iter < max_iter && improved && err_agg(end)>0

    improved = false;

    num_ext = size(external_edges, 1);
    to_check_max = ones(num_ext,1);
    improved1 = false;
    
    while ~improved1
        
        to_adjust = [];
        flag = 0;
        while ~flag
            [val_err_max,idx1] = max(external_edges(:,3));
            flag = to_check_max(idx1);
            if val_err_max <= 0 
                improved1 = true;
                break; 
            end 
            if ~flag 
                external_edges(idx1,3) = 0;
                to_adjust = [to_adjust;idx1];
            end
        end
        for i = to_adjust
            external_edges(i,3) = err(external_edges(i,1)) + err(external_edges(i,2));
        end
        if val_err_max <= 0 
            break; 
        end 
        
        i = external_edges(idx1,1); j = external_edges(idx1,2);
        
        found1 = false;
        tried = 0;  
        to_check = ones(num_ext,1);
        to_check(idx1) = 0;

        while ~found1
    
            flag = 0;
            while ~flag
                idx2 = randi([1,num_ext]);
                
                flag = to_check(idx2);
            end
            h = external_edges(idx2,1);
            k = external_edges(idx2,2);
            tried = tried + 1;
    
            if (i==h || j==k || i==k || j == h) 
                to_check(idx2) = 0;
            elseif A(i,h) || A(j,k) || A(i,k) || A(j,h)
                to_check(idx2) = 0;
            elseif length(unique([comm(i),comm(j),comm(h),comm(k)]))==4
                to_check(idx2) = 0;
            else
                cc1 = [ (comm(i) == comm(h)) ; (comm(j) == comm(k)) ];
                d1 = k_in_current([i,j,h,k]) + [cc1;cc1];
                e1 = k_in_desired([i,j,h,k]) - d1;
                e_1 = sum(e1);
                cc2 = [ (comm(i) == comm(k)) ; (comm(j) == comm(h)) ];
                d2 = k_in_current([i,j,k,h]) + [cc2;cc2];
                e2 = k_in_desired([i,j,k,h]) - d2;
                e_2 = sum(e2);
                [e_fin,idx] = min([e_1,e_2]);
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
            
            improved1 = true;
            count1 = count1 + 1;
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
            end
        end
        
        if ~improved1 
            to_check_max(idx1) = 0;
        end
    end
    
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
                improved2 = true;
                break; end
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
            elseif internal_edges(idx1,3) + internal_edges(idx1,3) >= -2
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
            count2 = count2 + 1;
            
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
            end
        end
        
        if ~improved2 
            to_check_min(idx1) = 0; 
        end
        
    end

    err_agg = [err_agg , mean(err)];
    var = (sum((err).^2))./N;
    var_agg = [var_agg , var];
    
    improved = improved1 || improved2;  
    
    iter = iter + 1;

end

A_rewired = A;

fprintf('[ iter          | %d \n  shuffle err>0 | %d \n  shuffle err<0 | %d \n  err           | %d ]\n', ...
    iter, count1, count2, mean(err));

end


