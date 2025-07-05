function [ci, H, bc] = fast_Stats(op, grp, nCommMax, edges)
% Calcola in **una sola passata** CI, entropia (normalizzata) e BC
% per tutte le comunità presenti in 'grp'.

    % ------- Pre-accumuli utili  --------------------------------------
    nBin   = numel(edges)-1;
    nNode  = accumarray(grp,1,[nCommMax 1]);             % nodi per comunità
    mu     = accumarray(grp,op,[nCommMax 1],@mean,NaN);   % media per comunità

    % 1) Concentration Index  (ε = 0.05 fisso qui, cambia se vuoi)
    epsCI  = 0.05;
    close  = abs(op - mu(grp)) <= epsCI;
    ci     = accumarray(grp,close,[nCommMax 1],@mean,NaN);

    % 2) Shannon entropy (25 bin)
    bin    = discretize(op, edges);
    hc     = accumarray([bin(:) grp(:)], 1, [nBin nCommMax]); % conteggi bin×comm
    p      = hc ./ nNode';                                    % frequenze
    p(p==0)= 1;                                               % evita log(0)
    H      = -sum( p .* log(p) ) ./ log(nBin);                % normalizzata

    % 3) Bimodality Coefficient
    %    Skewness & kurtosis via momenti centrali (vectorizzato)
    mu2    = accumarray(grp,(op - mu(grp)).^2,[nCommMax 1],@mean,NaN);
    mu3    = accumarray(grp,(op - mu(grp)).^3,[nCommMax 1],@mean,NaN);
    mu4    = accumarray(grp,(op - mu(grp)).^4,[nCommMax 1],@mean,NaN);
    g1     = mu3 ./ mu2.^(1.5);
    g2     = mu4 ./ mu2.^2 - 3;
    bc     = (g2 + 3) ./ (g1.^2 + 1);
end
