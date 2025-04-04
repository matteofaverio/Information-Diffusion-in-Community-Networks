function [W] = trustiness(A)

n = size(A,1);
W = zeros(n);

B = A | A';

deg = sum(B,2);
for i = 1:n
    neighbours = find(A(i,:));
    for j = neighbours
        p = deg(i)/(deg(i)+deg(j));
        q = deg(j)/(deg(i)+deg(j));
        x = rand;
        if p > x
            W(i,j) = p; W(j,i) = q;
        else
            W(i,j) = q; W(j,i) = p;
        end
        
    end
end

colSum = sum(W,1);
W = W./colSum;
W(isnan(W)) = 0;

end