function [finalOpinions, totalIterations, opinionHistory] = HK(A, W, initialOpinions, confidence)
% Optimized implementation of the HK model without changing its purpose.
%
% Inputs:
%   A - nxn adjacency matrix representing the network
%   W - nxn trustiness matrix
%   initialOpinions - nx1 vector of initial opinions
%   confidence - nx1 vector of confidence thresholds for each agent
%
% Outputs:
%   finalOpinions - nx1 vector of converged opinions
%   totalIterations - number of iterations until convergence
%   opinionHistory - matrix storing opinion evolution over time

%initialOpinions(agents) = 1;
n = length(initialOpinions);
tolerance = 1e-3;
maxIterations = 10000;
opinionChange = 1;
iteration = 0;
currentOpinions = initialOpinions;

% Precompute the neighbors for each agent from the adjacency matrix A.
neighborsCell = cell(n, 1);
for i = 1:n
    neighborsCell{i} = find(A(i, :));
end

% Preallocate the opinion history matrix.
opinionHistory = zeros(n, maxIterations + 1);
opinionHistory(:, 1) = currentOpinions;

while opinionChange > tolerance && iteration < maxIterations
    nextOpinions = zeros(n, 1);
    
    for i = 1:n

        neighbors = neighborsCell{i};
        % Identify influencers: neighbors whose opinions are within the confidence threshold.
        opinionDifferences = abs(currentOpinions(neighbors) - currentOpinions(i));
        withinConfidence = opinionDifferences < confidence(i);
        influencers = neighbors(withinConfidence);
        numInfluencers = length(influencers);
        
        if numInfluencers > 0
            % Normalize weights for the influencers.
            weights = W(influencers, i)./sum(W(influencers, i));
            % Calculate the weighted opinion
            weightedOpinion = sum(weights .* currentOpinions(influencers));
            % Combine self-opinion with the weighted average of influencers.
            nextOpinions(i) = (numInfluencers * weightedOpinion + currentOpinions(i)) / (numInfluencers + 1);
        else
            % If no influencer is within confidence, retain current opinion.
            nextOpinions(i) = currentOpinions(i);
        end
        %nextOpinions(agents) = 1;


    end
    
    % Calculate the total change in opinions and update the history.
    opinionChange = sum(abs(currentOpinions - nextOpinions));
    currentOpinions = nextOpinions;
    iteration = iteration + 1;
    opinionHistory(:, iteration + 1) = nextOpinions;
end

finalOpinions = currentOpinions;
totalIterations = iteration;
% Trim the opinion history to the number of iterations actually performed.
opinionHistory = opinionHistory(:, 1:(iteration + 1));
end
