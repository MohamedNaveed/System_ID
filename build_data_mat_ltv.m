function [V] = build_data_mat_ltv(U, y_mat, q, nu, nz, k, no_rollouts)


if k > q
    V = zeros((nz + nu)*q + nu, no_rollouts);

    V(1:nu,:) = U((k-1)*nu +1: k*nu,:);

    for i = 1:q

        V(nu + (i-1)*nz + 1: nu + i*nz, :) = y_mat((k - 1 - i)*nz + 1: (k-i)*nz, :);

        V(nu + q*nz + (i-1)*nu + 1:nu + q*nz + i*nu, :) = U((k - 1 - i)*nu +1: (k-i)*nu,:);
    end
    
else
    V = zeros((nz + nu)*(k-1) + nu, no_rollouts);

    V(1:nu,:) = U((k-1)*nu +1: k*nu,:);

    if k~=1
        for i = 1:k-1

            V(nu + (i-1)*nz + 1: nu + i*nz, :) = y_mat((k - 1 - i)*nz + 1: (k-i)*nz, :);

            V(nu + (k-1)*nz + (i-1)*nu + 1:nu + (k-1)*nz + i*nu, :) = U((k - 1 - i)*nu +1: (k-i)*nu,:);
        end
    end

end
end



