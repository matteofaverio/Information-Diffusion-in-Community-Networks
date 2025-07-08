function [H] = shannonEntropy(c, op, nbins, normalise)
%SHANNONENTROPY  Dispersione/frammentazione via entropia di Shannon.
%
%   [H, GROUPS] = shannonEntropy(c, op)                % 25 bin, normalizzata
%   [H, GROUPS] = shannonEntropy(c, op, nbins)         % normalizzata
%   [H, GROUPS] = shannonEntropy(c, op, nbins, normFlag)
%
%   INPUT
%     nbins      : numero di bin (default 20)
%     normalise  : true/false → se true divide per log(nbins)  (default true)
%
%   OUTPUT
%     H          : entropia (0 = tutto uguale, 1 = massima se normalizzata)

    if nargin < 3 || isempty(nbins),      nbins     = 25;   end
    if nargin < 4 || isempty(normalise),  normalise = true; end

    edges  = linspace(0, 1, nbins+1);
    groups = unique(c(:));
    H      = zeros(numel(groups),1);

    for k = 1:numel(groups)
        x       = op(c == groups(k));
        counts  = histcounts(x, edges);
        p       = counts / sum(counts);
        p       = p(p > 0);                 % evita log(0)
        h       = -sum(p .* log(p));
        if normalise
            h = h / log(nbins);
        end
        H(k) = h;
    end
end
