function [val,a,b] = find_couple_max(err,A)

    % A simmetrica: si considera solo la parte sopra-diagonale
    [i, j] = find(triu(A, 1));
    % Calcola la somma dei valori per ciascuna coppia connessa
    sums = err(i,2) + err(j,2);

    
    [val, idx] = max(sums);
    
    x = i(idx);
    y = j(idx);

    a = err(x,1);
    b = err(y,1);



end
