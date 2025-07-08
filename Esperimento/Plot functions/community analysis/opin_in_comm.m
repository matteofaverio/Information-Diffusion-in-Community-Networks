function [op_mean, op_mode, op_var] = opin_in_comm(c,op_end)
% op_mean -> vettore contenente in posizione i l'opinione media della
%            comunità i
% op_mode -> vettore contenente in posizione i la moda delle opinioni della
%            comunità i calcolata su 25 bin
% op_var -> vettore contenente in posizione i la varianza normalizzata
%           dell'opinione nella comunità i
% 
% 
% 
% 
% -
n = max(c);
op_mean = zeros(n,1);
op_var = zeros(n,1);
op_mode = zeros(n,1);
% definizione dei bin per il calcolo della moda
edges = linspace(0,1,26);                    

% Calcolo di media, moda e varianza dell'opinione in ciascuna comunità
for i = 1:n
    % trovo nodi della comunità i
    comm = (c == i);
    % trovo le opinioni finali della comunità i
    op_comm = op_end(comm);

    % calolco la media
    op_mean(i) = mean(op_comm);

    % calcolo la moda
    counts = histcounts(op_comm, edges);               % counts è 1×25
    [~, idxBin] = max(counts);                   % idxBin = indice del bin modale
    % per restituire il valore “modale” come centro del bin:
    binCenters = edges(1:end-1) + diff(edges)/2; % 1×25, centro di ogni bin
    op_mode(i) = binCenters(idxBin);

    % calcolo la varianza
    op_var(i) = var(op_comm)/0.25;
end


    
