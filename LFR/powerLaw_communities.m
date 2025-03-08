function [S,N, s_min, s_max] = powerLaw_communities(n,d_min,d_max, beta)
%powerLaw_communities restutuisce N numeri casuali distribuiti con legge di
%potenza troncata con esponente beta che sommano a n.
%
%   Input
%   n: somma del vettore
%   d_min: lower bound del minimo della distribuzione
%   d_max: upper bound del massimo della distribuzione
%   beta: esponente della distribuzione
%
%   Output
%   S:  vettore coi numeri casuali generati
%   N:  numero di campioni generati
%   s_min:  minimo della distribuzione
%   s_max:  massimo della distribuzione

% setting the minimum and the maximum size of the communities
s_min = d_min + randi([0 d_min]);
s_max = d_max + randi([0 round(d_max)]);
N = 0;
S = [];

tot = 0;
while tot < n - s_min
u = rand;
new = floor(( s_min^(1-beta) + u*(s_max^(1-beta) - s_min^(1-beta))).^(1/(1-beta)));
tot = tot + new;
N = N+1;
S(N) = new;
end

S(N) = n - sum(S(1:N-1));

end


