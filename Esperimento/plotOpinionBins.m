function plotOpinionBins(op, delta, step)
% plotOpinionBins Plots the evolution of node opinions grouped by initial bins.
%
%   plotOpinionBins(op) uses default delta = 0.002 and step = 1.
%   plotOpinionBins(op, delta) uses specified delta and default step = 1.
%   plotOpinionBins(op, delta, step) uses specified delta and step.
%
%   INPUTS:
%     op    - N×T matrix of opinions: op(i,t) is node i's opinion at time t (in [0,1]).
%     delta - (optional) bin width for initial opinions (default 0.002).
%     step  - (optional) sampling interval for time steps when plotting (default 1).
%
%   The function groups nodes into bins based on their initial opinion (column 1),
%   computes, for each bin, the average opinion at each time step, and plots
%   these series with colors proportional to the initial bin location.

% Handle default arguments
if nargin < 3 || isempty(step)
    step = 1;
end
if nargin < 2 || isempty(delta)
    delta = 0.002;
end

% Data dimensions
[N, T] = size(op);
% Time vector (sampled)	
tVec = 1:step:T;

% Bin assignment based on initial opinions
iniOpinions = op(:,1);
bin_idx = floor( iniOpinions./delta ) + 1;
K = max(bin_idx);

% Preallocate average-opinion matrix: K bins × T times
E = nan(K, T);
for b = 1:K
    nodes = (bin_idx==b);
    if ~any(nodes)
        continue;
    end
    E(b,:) = mean(op(nodes,:), 1);
end

% Plotting
figure;
hold on;
cm = jet(K);  % colormap
for b = 1:K
    y = E(b, tVec);
    h = plot(tVec, y, 'LineWidth', 1.2, 'Color', cm(b,:));
    % Add transparency if MATLAB supports it
    if verLessThan('matlab', '9.5') == 0
        h.Color(4) = 0.3;  % alpha channel
    end
end

xlabel('Passi temporali');
ylabel('Opinione media del bin');
title('Evoluzione delle opinioni raggruppate per bin iniziale');
grid on;

% Colorbar
c = colorbar('Ticks', [0, 0.25, 0.5, 0.75, 1], ...
             'TickLabels', {'0', '0.25', '0.5', '0.75', '1'});
c.Label.String = 'Opinione iniziale';
colormap(cm);

hold off;
end
