function [op_mean,op_var] = opin_in_comm(c,op_end)

n = max(c);
op_mean = zeros(n,1);
op_var = zeros(n,1);

for i = 1:n
    comm = find(c==i);
    op_comm = op_end(comm);
    op_mean(i) = mean(op_comm);
    op_var(i) = var(op_comm);
end


