function [A,weighted_err] = rewiring_LFR(A,dd,c,mu,maxit,rs)
%%
n = size(A,1);
it = 0;

% Calcola current_int per ogni nodo e l'errore
current_int = zeros(n,1);
for i = 1:n
    for j = 1:n
        if (c(i) == c(j) && A(i,j) == 1)
            current_int(i) = current_int(i) + 1;
        end
    end
end

dd_real = sum(A,2);
abs_err = mu - current_int./dd_real;
weighted_err = abs_err.*dd_real;
diff = dd_real - dd;


%%
x=1;
count = 0;
to_inspect1 = ones(n,1);


while( it<maxit && x>0 )

    flag=0;

    while( flag==0 )
        % trovo l'indice del nodo con errore massimo
        [x,i] = max(weighted_err);
        if ( x<=0 )
            flag = 1; % se non esistono più nodi con errore positivo da sistemare, esci dal ciclo
        else
            % creo V_i, il set dei nodi adiacenti a i e di altre comunità e 
            % W_i, il set dei nodi non adiacenti a i e della stessa comunità
            V_i = intersect( find(A(i,:)==1) , find(c~=c(i)) );
            A(i,i) = 1; % evito che i sia inserito in W_i
            W_i = intersect( find(A(i,:)==0) , find(c==c(i)) );
            A(i,i) = 0; % tolgo il self loop 
            if (isempty(V_i) || isempty(W_i))% se almeno uno dei due set è vuoto, ignoro i
                weighted_err(i) = 0;
                to_inspect1(i) = 0;
            else % cerco in V_i il nodo da staccare e in W_i quello da attaccare
                [~,idx1] = max(weighted_err(V_i));
                j = V_i(idx1);
                [~,idx2] = max(weighted_err(W_i));
                h = W_i(idx2);
                flag = 1; % ho trovato i due candidati, posso uscire dal ciclo
            end
        end
    end

    if( x>0 ) % eseguo solo se ho trovato i tre candidati
        
        % stacco i e j, attacco i e h
        A(i,j) = 0; A(j,i) = 0;
        A(i,h) = 1; A(h,i) = 1;

        % aggiorno current_int, dd_real e diff
        current_int(i) = current_int(i) + 1; current_int(h) = current_int(h) + 1;       
        dd_real(j) = dd_real(j) - 1; dd_real(h) = dd_real(h) + 1;
        diff(j) = diff(j) - 1; diff(h) = diff(h) + 1;

        if(rs == 1)
            flag1=0;
            while( flag1==0 )
                [diff_p,p] = max(diff);
                if ( diff_p<=0 )
                    flag1 = 1;
                else
                    V_p = find(A(p,:)==1);
                    if (isempty(V_p))
                        diff(p) = 0;
                    else
                        R = intersect( max(diff(V_p)), diff~=0);
                        if isempty(R)
                            diff(p) = 0;
                        else
                            [~,idx_r] = max(diff(V_p));
                            r = V_p(idx_r);
                            A(p,r) = 0; A(r,p) = 0; 
                            dd_real(p) = dd_real(p) - 1; dd_real(r) = dd_real(r) - 1;
                            diff(p) = diff(p) - 1; diff(r) = diff(r) - 1;
                            flag1 = 1;
                        end
                    end
                end
            end 
    
            flag2=0;
            while( flag2==0 )
                [diff_s,s] = min(diff);
                if ( diff_s>=0 )
                    flag2 = 1;
                else
                    V_s = find(A(s,:)==0);
                    if (isempty(V_s))
                        diff(s) = 0;
                    else
                        T = intersect( min(diff(V_s)), diff~=0);
                        if isempty(T)
                            diff(s) = 0;
                        else
                            [~,idx_t] = max(diff(V_s));
                            t = V_s(idx_t);
                            A(s,t) = 1; A(t,s) = 1; 
                            dd_real(s) = dd_real(s) + 1; dd_real(t) = dd_real(t) + 1;
                            diff(s) = diff(s) - 1; diff(t) = diff(t) - 1;
                            current_int(s) = current_int(s) + 1; 
                            current_int(t) = current_int(t) + 1;
                            flag2 = 2;
                        end
                    end
                end
            end
        end 
    
    diff = dd_real - dd;
    weighted_err = (mu*dd_real - current_int).*to_inspect1; % ha senso inspect?

    end
    
    
    it = it + 1;
    
    if mean(mu*dd - current_int) < 0 % se sono vicino al mu desiderato, mi fermo
            x = 0
    end



end


it