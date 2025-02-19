function [finalState, totalSteps, activationHistory] = LTM(adjacencyMatrix, initialState, thresholdVector)

%
% Linear Threshold Model (LTM) Simulation Until Convergence
%
% Inputs:
%   - adjacencyMatrix: nxn adjacency matrix representing the network (A)
%   - initialState: nx1 binary vector representing active (1) and inactive (0) nodes (H_t0)
%   - thresholdVector: nx1 vector specifying the activation threshold for each node
%
% Outputs:
%   - finalState: nx1 binary vector representing the state of nodes after diffusion (H_final)
%   - totalSteps: Number of steps taken until convergence
%   - activationHistory: Matrix storing the state of nodes at each step
%


numNodes = size(adjacencyMatrix, 1); % Number of nodes
currentState = initialState; % Initialize the state
activationHistory = currentState; % Store states over time
totalSteps = 0; % Count the number of diffusion steps

while true
    nextState = currentState; % Initialize next state as the current state
    newlyActivated = 0; % Track new activations in this step
    
    for node = 1:numNodes
        if currentState(node) == 0 % Only check inactive nodes
            influenceSum = sum(adjacencyMatrix(:, node) .* currentState); % Compute influence
            
            if influenceSum >= thresholdVector(node) % Check activation condition
                nextState(node) = 1; % Activate node
                newlyActivated = newlyActivated + 1;
            end
        end
    end
    
    totalSteps = totalSteps + 1; % Increment step count
    activationHistory(:, totalSteps + 1) = nextState; % Store current state
    
    if newlyActivated == 0 % Stop if no new activations occur
        break;
    end
    
    currentState = nextState; % Update state for next iteration
end

finalState = currentState; % Final state after diffusion ends

end


