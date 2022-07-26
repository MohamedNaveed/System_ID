function [V] = build_V(u, y, p)

m = size(u,1);
l = size(u,2);
q = size(y,1);

V = zeros((m+q)*p + m, l);

V(1:m,:) = u;

for i = 1:p
   
    V(m + (i-1)*(m+q) + 1: m + i*(m+q),:) = [zeros(m,i) u(:,1:l-i);...
                                            zeros(q,i) y(:,1:l-i)];
    
end



end

