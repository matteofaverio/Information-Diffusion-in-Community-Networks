function [CI] = concentrationIndex_opt(c, op, delta)
% Versione verificata e robusta di concentrationIndex.
if nargin < 3,  delta = 0.05;  end

% VALIDAZIONE E FILTRAGGIO
c = c(:);
op = op(:);
if numel(c) ~= numel(op)
    error('concentrationIndex_opt:InputSizeMismatch', ...
          'I vettori di input c (comunità) e op (opinioni) devono avere lo stesso numero di elementi.');
end
valid_idx = ~isnan(op);
op_valid = op(valid_idx);
c_valid  = c(valid_idx);

if isempty(c_valid)
    CI = [];
    return;
end

% Logica di calcolo (corretta)
[group_idx] = findgroups(c_valid);
fcn = @(x) mean(abs(x - mean(x)) <= delta);
CI = splitapply(fcn, op_valid, group_idx);
end