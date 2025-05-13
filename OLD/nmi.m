function nmi_value = nmi(c1, c2)
% nmi calcola la Normalized Mutual Information (NMI) tra due partizioni.
%
% USO:
%   nmi_value = nmi(nodes, c1, c2)
%
% INPUT:
%   nodes - vettore contenente l'elenco dei nodi.
%   c1    - vettore che rappresenta la partizione dei nodi (prima partizione),
%           dove c1(i) = j se il nodo i appartiene al cluster j.
%   c2    - vettore che rappresenta la partizione dei nodi (seconda partizione),
%           dove c2(i) = j se il nodo i appartiene al cluster j.
%
% OUTPUT:
%   nmi_value - valore della normalized mutual information tra le due partizioni.
%
% La NMI è calcolata come:
%   NMI = I(c1, c2) / sqrt(H(c1)*H(c2))
% dove I(c1, c2) è la mutua informazione e H(c1), H(c2) sono le entropie
% delle partizioni.

N = length(c1);

% Assicurarsi che i vettori di partizione siano colonne
c1 = c1(:);
c2 = c2(:);

% Trova i cluster unici per ciascuna partizione
clusters1 = unique(c1);
clusters2 = unique(c2);

% Calcolo dell'entropia per la prima partizione
H1 = 0;
for i = 1:length(clusters1)
    p_i = sum(c1 == clusters1(i)) / N;
    if p_i > 0
        H1 = H1 - p_i * log(p_i);
    end
end

% Calcolo dell'entropia per la seconda partizione
H2 = 0;
for j = 1:length(clusters2)
    p_j = sum(c2 == clusters2(j)) / N;
    if p_j > 0
        H2 = H2 - p_j * log(p_j);
    end
end

% Calcolo della mutua informazione (MI)
MI = 0;
for i = 1:length(clusters1)
    for j = 1:length(clusters2)
        % Numero di nodi in comune tra il cluster i in c1 e il cluster j in c2
        n_ij = sum((c1 == clusters1(i)) & (c2 == clusters2(j)));
        p_ij = n_ij / N;
        if p_ij > 0
            p_i = sum(c1 == clusters1(i)) / N;
            p_j = sum(c2 == clusters2(j)) / N;
            MI = MI + p_ij * log(p_ij / (p_i * p_j));
        end
    end
end

% Calcolo della normalized mutual information
if H1 == 0 || H2 == 0
    if H1 == 0 && H2 == 0
        nmi_value = 1; % entrambe le partizioni sono degenere e identiche
    else
        nmi_value = 0;
    end
else
    nmi_value = MI / sqrt(H1 * H2);
end

end


