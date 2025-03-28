function [W] = trustiness(A)
n = size(A,1);
W = zeros(n);

B = A | A';

deg = sum(B,2);
for i = 1:n
    for j = 1:n
        if i~=j
            total = deg(i)+deg(j);
            w = (total*rand - deg(i))/total;
                if w > 0
                    W(i,j) = w;
                else
                    W(j,i) = -w;
                end
        end
        
    end
end
W = A.*W;
colSum = sum(W,1);
W = W./colSum;
W(isnan(W)) = 0;

end