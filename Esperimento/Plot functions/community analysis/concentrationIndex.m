function [CI] = concentrationIndex(c, op, delta)
%CONCENTRATIONINDEX  Consenso locale: quota di nodi entro ±epsilon dalla media.
%
%   [CI, GROUPS] = concentrationIndex(c, op)
%   [CI, GROUPS] = concentrationIndex(c, op, epsilon)
%
%   INPUT
%     c       : vettore comunità (1×N o N×1)
%     op      : vettore opinioni (stesso size di c)
%     delta : raggio (default 0.05)
%
%   OUTPUT
%     CI      : array (Nc×1) con il Concentration Index di ogni comunità
%     GROUPS  : valori distinti di c nell’ordine restituito

    if nargin < 3,  delta = 0.05;  end

    groups = unique(c(:));
    CI     = zeros(numel(groups),1);

    for k = 1:numel(groups)
        idx   = (c == groups(k));
        x     = op(idx);
        mu    = mean(x);
        CI(k) = mean(abs(x - mu) <= delta);
    end
end
