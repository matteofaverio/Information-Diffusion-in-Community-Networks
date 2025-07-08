function [BC] = bimodalityCoefficient_opt(c, op)
% Versione verificata e robusta di bimodalityCoefficient.

% VALIDAZIONE E FILTRAGGIO
c = c(:);
op = op(:);
if numel(c) ~= numel(op)
    error('bimodalityCoefficient_opt:InputSizeMismatch', ...
          'I vettori di input c (comunità) e op (opinioni) devono avere lo stesso numero di elementi.');
end
valid_idx = ~isnan(op);
op_valid = op(valid_idx);
c_valid  = c(valid_idx);

if isempty(c_valid)
    BC = [];
    return;
end

% Logica di calcolo (corretta)
[group_idx] = findgroups(c_valid);
fcn = @(x) kurtosis(x, 1) / (skewness(x, 1)^2 + 1);
BC = splitapply(fcn, op_valid, group_idx);
end