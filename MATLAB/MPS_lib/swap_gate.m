function U = swap_gate(d)
%swap_gate gives a unitary matrix implementing swap gate on two qudits
%   
%   U = swap_gate(d)
%       U = \sum_{ij} |ij><ji|
%
%   Note U is a d^2 x d^2 dimensional unitary matrix

zs = eye(d);

U = 0;
for ind1 = 1:d
    z1 = zs(:,ind1);
    
    U = U + diag(kron(z1, z1));
    
    for ind2 = ind1+1:d
        z2 = zs(:, ind2);
        U = U + kron(z1, z2) * kron(z2, z1)' + kron(z2, z1) * kron(z1, z2)';
    end
end

end