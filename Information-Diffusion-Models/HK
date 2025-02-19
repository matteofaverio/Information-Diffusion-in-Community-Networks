function [finalOpinions, totalIterations, opinionHistory] = HK(adjacencyMatrix, initialOpinions, confidenceLevels)

%
% Inputs:
%   - adjacencyMatrix: nxn adjacency matrix representing the network (A)
%   - initialOpinions: nx1 vector representing initial opinions (x_0)
%   - confidenceLevels: nx1 vector specifying the confidence threshold for each agent
%
% Outputs:
%   - finalOpinions: nx1 vector of converged opinions
%   - totalIterations: Number of iterations until convergence
%   - opinionHistory: Matrix storing opinion evolution over time
%
% Usage:
%   [finalOpinions, totalIterations, opinionHistory] = HK_model(A, x_0, confidence)

n = length(initialOpinions);
opinionChange = 1;  % Track opinion changes
tolerance = 1e-3;   % Convergence threshold
maxIterations = 100; % Maximum allowed iterations
iteration = 0;      % Iteration counter

currentOpinions = initialOpinions;
opinionHistory = currentOpinions; % Store opinions at each step

while (opinionChange > tolerance && iteration < maxIterations)
    
    % Compute the neighborhood matrix based on confidence levels and adjacency
    I = (abs(currentOpinions - currentOpinions') < confidenceLevels) & adjacencyMatrix;
    
    % Compute next opinions 
    neighborCounts = sum(I, 2); % Number of neighbors per agent, including itself
    nextOpinions = currentOpinions; % Initialize next state
    validAgents = (neighborCounts > 1); % Agents that have neighbors
    nextOpinions(validAgents) = sum(I(validAgents, :) .* currentOpinions', 2) ./ neighborCounts(validAgents);
    
    % Store opinion history
    opinionHistory = [opinionHistory, nextOpinions];

    % Compute total opinion change
    opinionChange = sum(abs(currentOpinions - nextOpinions));

    % Update opinions and iteration count
    currentOpinions = nextOpinions;
    iteration = iteration + 1;
end

finalOpinions = currentOpinions;
totalIterations = iteration;

end
