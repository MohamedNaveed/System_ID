function [P] = calculate_P(h_o, k, q, nz)


P = zeros(q*nz,nz);

for j = 1:q

    P((j-1)*nz +1 : j*nz,:) = h_o((k+j-1)*nz + 1: (k+j)*nz,(k-1)*nz + 1:(k)*nz);  
end
    
end

