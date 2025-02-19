%% Linear Threshold Model (LTM) for Information Diffusion
% This function simulates one step of the Linear Threshold Model (LTM) in a network.
% 
% Inputs:
%   - adjacencyMatrix: nxn adjacency matrix representing the network (A)
%   - initialState: nx1 binary vector representing active (1) and inactive (0) nodes 
% 
% Outputs:
%   - nextState: nx1 binary vector representing the state of nodes after diffusion 
%   - influencedNodes: Number of newly activated nodes in this step (k)
% 

function [nextState, influencedNodes] = LTM(adjacencyMatrix, initialState)

numNodes = size(adjacencyMatrix, 1); % Number of nodes in the network
influencedNodes = 0; % Counter for newly influenced nodes
nextState = initialState; % Initialize next state as the current state

for node = 1:numNodes
    if initialState(node) == 0 % Consider only inactive nodes
        influenceSum = 0; % Sum of incoming influences
        
        for neighbor = 1:numNodes
            influenceSum = influenceSum + adjacencyMatrix(neighbor, node) * initialState(neighbor); 
            % Sum influences from already active neighbors
        end
        
        if influenceSum >= adjacencyMatrix(node, node) % Check activation threshold
            influencedNodes = influencedNodes + 1; % Increase the count of newly activated nodes
            nextState(node) = 1; % Activate the node
        end
    end
end

end

