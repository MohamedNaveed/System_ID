function [P] = calculate_P(h_o, k, q, nz)


P = zeros(q*nz,nz);

for j = 1:q
    
    rows = (k+j-1)*nz + 1: (k+j)*nz;
    cols = (j-1)*nz + 1:(j)*nz;
    P((j-1)*nz +1 : j*nz,:) = h_o(rows, cols);  
end
    
end

