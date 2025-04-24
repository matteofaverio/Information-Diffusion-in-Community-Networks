function avgFraction = mu(A,c)

% maschera booleana, M(i,j) == 1 se i e j fanno parte della stessa comunità
M = (c == c');

% applico la maschera ad A e sommo sulle righe: trovo i collegamenti
% interni di ogni nodo
sameCommCounts = sum(A .* M, 2);

% trovo i collegamenti totali dei nodi 
degrees = sum(A, 2);
fractions = sameCommCounts ./ degrees;

% tolgo i nan
fractions(degrees == 0) = 0;
avgFraction = mean(fractions);

end

