function [opinions,t] = opinionProp_IC(W,v0)

%opinionProp_IC simulates opinion propagation in a network using
%Independent Cascade model
%   [opinions, t] = opinionProp_IC(W,v0)
%   Inputs:
%   W: trustiness matrix
%   v0: initial seed
%
%   Outputs:
%   opinions: propagated opinions at every time step from v0 to t
%   t: number of time steps

% Creazione del vettore contenente le opinioni nel tempo
opinions = v0;

% Vettore contenente i nodi appena attivati:
newly_activated = v0;

% inizializzo il tempo
t = 1;

opinions_updated = v0;

while sum(newly_activated) > 0

    for i = 1:length(v0)
        if opinions_updated(i) == 0
            opinions_updated(i) = opinions_updated(i) + binornd(1,W(:,i))'*newly_activated;
            if opinions_updated(i) ~= 0
                opinions_updated(i) = opinions_updated(i)/opinions_updated(i); 
            end
        end
    end

% aggiorno il vettore con i nuovi agenti attivi
newly_activated = opinions_updated - opinions(:,t);

% aggiorno l'evoluzione degli agenti attivi
if sum(newly_activated) > 0
    opinions = [opinions, opinions_updated];
end
% faccio avanzare il tempo
t = t+1;

end
end

