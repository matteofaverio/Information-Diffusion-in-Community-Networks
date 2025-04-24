function [dd,d_max] = powerLaw_degree(n,gamma,d_min,d)
%powerLaw_degree restituisce n numeri casuali con distribuzione in legge di
%potenza troncata di esponente gamma e media d. gli estremi della
%distribuzione sono dati da d_min, passato come input, e d_max, calcolato
%dalla funzione 
%
%   Input:
%   n: numero di generazioni
%   gamma: esponente della legge di potenza
%   d_min: lower bound
%   d: media del campione
%
%   Output
%   dd: vettore coi campioni casuali
%   d_max: upper bound


% la funzione che definisce il valore atteso della distribuzione presenta
% una singolarità eliminabile in gamma = 2, singolaritàa non gestibile dal
% calcolatore
if gamma == 2
    gamma = gamma + 1e-12;
end

% calcolo di d_max con fissato

num = @(x) (x.^(2-gamma)-d_min.^(2-gamma))/(2-gamma);
den = @(x) (x.^(1-gamma)-d_min.^(1-gamma))/(1-gamma);
dmean = @(x) num(x)./den(x);

% funzione a cui applicare il metodo di bisezione
f = @(x) dmean(x) - d;

% estremi dell'intervallo su cui applicare la bisezione
lower = d_min + 1e-9;
upper = 100.0;

% sposto l'upper bound finchè non ottengo un intervallo adeguato su cui
% applicare il metodo di bisezione.
while f(lower) * f(upper) > 0
    upper = upper * 2;
    % mostro un errore se l'intervallo diventa troppo ampio.
    if upper > 1e9
        error(['Non riesco a trovare un intervallo adeguato per d_max,' ...
            ' Prova ad altri parametri o a ridurre la media desiderata.']);
    end
end
% uso il metodo di bisezione per trovare d_max;
d_max = fzero(f,[lower, upper]);
d_max = round(d_max);

% Generazione dei gradi dei singoli nodi
u = rand(n,1);
dd = floor(( d_min^(1-gamma) + u*(d_max^(1-gamma) - d_min^(1-gamma))).^(1/(1-gamma)));

end

