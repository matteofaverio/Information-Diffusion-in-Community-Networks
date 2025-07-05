function [BC] = bimodalityCoefficient(c, op)
%BIMODALITYCOEFFICIENT  Indice di bimodalità basato su skewness e kurtosi.
%
%   [BC, GROUPS] = bimodalityCoefficient(c, op)
%
%   Formula usata:
%       BC = (gamma2 + 3) / (gamma1^2 + 1)
%   dove  gamma1 = skewness,  gamma2 = kurtosi in eccesso.

    groups = unique(c(:));
    BC     = zeros(numel(groups),1);

    for k = 1:numel(groups)
        x   = op(c == groups(k));
        g1  = skewness(x,0);        % bias-corrected
        g2  = kurtosis(x,0) - 3;    % eccesso
        BC(k) = (g1^2 + 1) / (g2 + 3);  
    end
end
