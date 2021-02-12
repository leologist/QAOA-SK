function QMPS = getQb_as_MPS(p, betas)
%getQb_as_MPS returns Q_b as a rank-1 MPS
%
%    QMPS = getQb_as_MPS(p, betas)

d = 4;

bs = [1     1;...
      1    -1;...
     -1     1;...
     -1    -1];
 
b_pj = bs(:, 1);
b_mj = bs(:, 2);

L_Q = cell(p-1,1);
G_Q = cell(p,1);


for j = 1:p-1
    L_Q{j} = 1;
end

for j = 1:p
    G_Q{j} = zeros(1,d);
    G_Q{j}(1,:) = cos(betas(j)).^(1+ (b_pj + b_mj)/2)  ...
                  .* sin(betas(j)).^(1- (b_pj + b_mj)/2) ...
                   .* (1i).^((b_mj - b_pj)/2);
   if j >= 2 && j <= p-1
       G_Q{j} = reshape(G_Q{j}, 1,1,d);
   end
end

QMPS.G = G_Q;
QMPS.L = L_Q;