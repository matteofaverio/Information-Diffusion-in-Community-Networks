function [H_t1,k] = LTM(A , H_t0)

n = size(A);
k = []; % numero nodi influenzati
H_t1 = H_t0;


for i = 1:n
    if H_t0(i) == 0 % considero i nodi non influenzati
        s = 0; % somma degli attacchi dai nodi adiacenti
        for j = 1:n
            s = s + A(j,i)*H_t0(j); % sommo gli attacchi da tutti i nodi,
                                    % considerando solo quelli già attivi
        end
        if s >= A(i,i)
            k = k+1; % un nuovo nodo influenzato
            flag = 1;
            H_t1(i) = 1;
        end
    end
end
