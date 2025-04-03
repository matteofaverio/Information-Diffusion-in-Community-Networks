function [finalOpinions, it, opinionHistory] = HK(A, W, initialOpinions, confidence)

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
%

n = length(initialOpinions);
opinionChange = 1;  % Track opinion changes
toll = 1e-3;   % Convergence threshold
maxit = 100; % Maximum allowed iterations
it = 0;      % Iteration counter

currentOpinions = initialOpinions;
opinionHistory = currentOpinions; % Store opinions at each step

while (opinionChange > toll && it < maxit)
    
    nextOpinions = zeros(n,1);

    for i = 1:n
        flag = 1
        neighbours = find(A(i,:));
        close_enough = find( abs(currentOpinions(neighbours)-currentOpinions(i)) < confidence(i));
        influencers = neighbours(close_enough);
        I = length(influencers);

        nextOpinions(i) = (W(influencers,i).*currentOpinions(influencers));
        nextOpinions(i) = (I*nextOpinions(i) + currentOpinions(i))/(I+1);

    end

    opinionChange = sum(abs(currentOpinions - nextOpinions));
    currentOpinions = nextOpinions;
    opinionHistory = [opinionHistory, nextOpinions];
    it = it + 1;
    
end

finalOpinions = currentOpinions;

end

