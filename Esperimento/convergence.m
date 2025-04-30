function n = convergence(X,toll)
%CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here

% last = X(:,end);
% res = max(abs(X-last));
% n = find(res>1e-03);
% n = n(end);

diff = max(abs(X(:,2:end)-X(:,1:end-1)));
n = find(diff>toll);
n = n(end);

end

