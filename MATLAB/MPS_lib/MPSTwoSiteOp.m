function [MPS, err_trunc] = MPSTwoSiteOp(MPS, Op, k, renormalize)
%MPSTwoSiteOp apply operator Op on sites k and k+1
%
%   [MPS, err_trunc] = MPSTwoSiteOp(MPS, Op, k, renormalize)
%
%   Input:
%     MPS = struct containing MPS in canonical form
%           required to have two fields:
%             MPS.G = a cell array containing N Gamma_{ab}^s tensors
%             MPS.L = a cell array containing N-1 diagonal matrices of
%                     Schmidt values
%     Op = d^2 x d^2 matrix to be applied to two sites (k, k+1)
%     k = location
%     renormalize = 1 if want output to represent normalized state
%                 Op*|psi>/||Op*|psi>||, 0 if just want Op*|psi>
%
%   Output:
%     MPS = state after applying two-site operation
%          (inherits all other struct fields beyond (G, L)of input MPS)
%     err_trunc = the error from truncating Schmidt values
%
%   This function assumes d, maxD, seed are already defined globally
%
%   Based on Vidal (2003), arXiv:quant-ph/0301063
%   Author: Leo Zhou
%   Revised: Feb 12, 2021

N_sites = length(MPS.G);
% check input
if ~(k>=1 && k < N_sites && length(MPS.L) == N_sites - 1)
    error('Invalid input MPS');
end

global d maxD seed

Gs = MPS.G;
Ls = MPS.L;
G_out = Gs;
L_out = Ls;
err_trunc = 0;

min_schmidt_val = 1e-6;

% optional parameters

if k == 1
    DL = 1;
    DM = length(Ls{1});
    if N_sites == 2
        DR = 1;
    else
        DR = length(Ls{2});
    end

    G1 = Gs{1}.'*Ls{1};
    if N_sites ~= 2
        G2 = reshape(permute(Gs{2}, [1,3,2]), [DM*d, DR])*Ls{2};
    else
        G2 = reshape(Gs{2}, [DM*d, DR]);
    end
    G2 = reshape(full(G2), [DM, d*DR]);

    tensor = reshape(G1*G2, [d, d, DR]); % i, j, a

    tensor = reshape(permute(tensor, [2,1,3]), [d^2, DR]); % (ji, a)
    tensor = Op * tensor; %Op = (jpip), (ji)
    tensor = reshape(permute(reshape(tensor, [d d, DR]), [2,1,3]), [d, d*DR]); % i, (ja)
        
%     tensor = reshape(G1*G2, [d^2, DR]);
%     tensor = Op * tensor;
%     tensor = reshape(tensor, [d d*DR]); % i, (ja)
    
    if renormalize
        tensor = renorm(tensor);
    end

    [U, S, V] = schimdt_decomp(tensor); % U = dxd, S = d x dDR, V = dDR x dDR
    S = diag(S);
    DM = nnz(S>min_schmidt_val);
%     L_out{1} = sparse(diag(S(1:DM)));
    L_out{1} = spdiags(S(1:DM), 0, DM, DM);
    
    %U_{i,gamma} S_{gamma,gamma} V_{gamma, ja}
    
    % G1 needs to be a DMxd tensor
    G_out{1} = transpose(U(:, 1:DM));
    
    % G2 needs to be DMxDRxd tensor
                                        %j  a  gamma
    V = permute(conj(reshape(V(:, 1:DM), [d, DR, DM])), [2, 1, 3]);
    V = reshape(V, [DR, d*DM]); % a, (jgamma)
    if N_sites > 2
        V = effinv(Ls{2})*V;
    end
    
                                  %a  j  gamma
    G_out{2} = permute(reshape(V, [DR, d, DM]), [3, 1, 2]);
    if N_sites == 2
        G_out{2} = reshape(G_out{2}, [DM, d]);
    end
    
elseif k == N_sites-1
    DL = length(Ls{k-1});
    DM = length(Ls{k});
    DR = 1;
%     tensor = zeros(DL, d, d);
%     for phys_ind1 = 1:d
%         for phys_ind2 = 1:d
% %             size(Lambdas{k-1})
%             tensor(:, phys_ind1, phys_ind2) = ...
%                 Ls{k-1} * Gs{k}(:,:, phys_ind1) * Ls{k} * Gs{k+1}(:, phys_ind2);
%                 %  DLxDL            DLxDM                       DMxDM            DMx1
%         end
%     end
%     tensor = reshape(permute(tensor, [2 3 1]), [d^2, DL]);

    G1 = Ls{k-1}*reshape(Gs{k}, [DL, DM*d]);
    G1 = reshape(permute(reshape(full(G1),[DL, DM, d]), [1,3,2]), [DL*d, DM]);
    G2 = Ls{k}*Gs{k+1};
    
    
    tensor = reshape(G1*G2, [DL, d, d]); % a, i, j
    tensor = reshape(permute(tensor, [3,2,1]), [d^2, DL]); % (ji, a)

    tensor = Op * tensor; %Op = (jpip), (ji)
                                             %(j, i, a)   
    tensor = reshape(permute(reshape(tensor, [d, d, DL]), [3,2,1]), [DL*d, d]); % (ai),j
    
%     tensor = reshape(G1*G2, [DL, d^2]).';
%     tensor = Op * tensor;
%     tensor = reshape(tensor, [d d DL]); % i j a
%     tensor = permute(tensor, [3 1 2]); % a i j
%     tensor = reshape(tensor, [DL*d, d]); %(ai), j
    
    if renormalize
        tensor = renorm(tensor);
    end
    
    [U, S, V] = schimdt_decomp(tensor); %U = DLd x DLd, S = DLd x d, V = dxd
    S = diag(S);
    DM = nnz(S>min_schmidt_val);
%     L_out{k} = sparse(diag(S(1:DM)));
    L_out{k} = spdiags(S(1:DM), 0, DM, DM);
    
                           %a  (igamma)
    U = reshape(U(:, 1:DM), [DL, d*DM]); % truncate
    U = effinv(Ls{k-1})*U;
                                  %a  i gamma
    G_out{k} = permute(reshape(U, [DL, d, DM]), [1, 3, 2]); % DxDxd tensor
    
    G_out{k+1} = transpose(conj(V(:, 1:DM))); % DMxd tensor   
    
else
    DL = length(Ls{k-1});
    DM = length(Ls{k});
    DR = length(Ls{k+1});
%     tensor = zeros(DL, DR, d, d);
%                   %a   b   i  j
%     for phys_ind1 = 1:d
%         for phys_ind2 = 1:d
%             tensor(:, :, phys_ind1, phys_ind2) = ...
%                 Ls{k-1} * Gs{k}(:, :, phys_ind1) * Ls{k} * Gs{k+1}(:, :, phys_ind2) * Ls{k+1};
%                %   DLxDL            DLxDM          DMxDM           DMxDR               DRxDR
%         end
%     end
%     tensor = reshape(permute(tensor, [3 4 1 2]), [d^2, DL*DR]);
    G1 = reshape(permute(Gs{k}, [1,3,2]), [DL, d*DM]);
    G1 = reshape(full(Ls{k-1}*G1), [DL*d, DM]);
    G2 = reshape(permute(Gs{k+1}, [1,3,2]), [DM*d, DR]);
    G2 = reshape(full(G2*Ls{k+1}), [DM, d*DR]);

    
                                                  % a, i, j, b     (j, i,  a, b
    tensor = reshape(permute(reshape(G1*Ls{k}*G2, [DL, d, d, DR]), [3, 2, 1, 4]), [d^2, DL*DR]);
    
    tensor = Op * tensor; %(jpip,ji)
    tensor = reshape(tensor, [d d DL DR]); % j i a b
    tensor = permute(tensor, [3 2 4 1]); % a i b j
    tensor = reshape(tensor, [DL*d, DR*d]); % (ai), (bj)
    
%     tensor = reshape(permute(reshape(G1*Ls{k}*G2, [DL, d, d, DR]), [2, 3, 1, 4]), [d^2, DL*DR]);
%     tensor = Op * tensor;
%     tensor = reshape(tensor, [d d DL DR]); % i j a b
%     tensor = permute(tensor, [3 1 4 2]); % a i b j
%     tensor = reshape(tensor, [DL*d, DR*d]); % (ai), (bj)
    
    if renormalize
        tensor = renorm(tensor);
    end
    
    [U, S, V] = schimdt_decomp(tensor);
    S = diag(S);
    if ~issorted(wrev(S))
        error('Singular values not sorted?!');
    end
    DM = nnz(S>min_schmidt_val);
    if DM > maxD %truncate
        DM = maxD;
    end
    S_norm = norm(S);
    S_trunc = S(1:DM);
    err_trunc = abs(S_norm^2-norm(S_trunc)^2)/S_norm^2;
    if norm(S_trunc) > 0
        S_trunc = S_trunc/norm(S_trunc)*S_norm;
    end
%     L_out{k} = sparse(diag(S_trunc));
    L_out{k} = spdiags(S_trunc, 0, DM, DM);
    
                           %a  (igamma)
    U = reshape(U(:, 1:DM), [DL, d*DM]); % truncate
    U = effinv(Ls{k-1})*U;
    U = reshape(U, [DL, d, DM]);
    G_out{k} = permute(U, [1, 3, 2]); %a gamma i
    
                                %b  (jgamma)
    V = conj(reshape(V(:, 1:DM), [DR, d*DM])); %truncate
    V = effinv(Ls{k+1})*V;
    V = reshape(V, [DR, d, DM]);
    G_out{k+1} = permute(V, [3, 1, 2]); %gamma b j
end

MPS.G = G_out;
MPS.L = L_out;


    function [U, S, V] = schimdt_decomp(tensor)
        try
            [U, S, V] = svdecon(tensor); %schmidt decomposition 
        catch ME
            fprintf(2, 'SVD error: seed = %i,   site ind = %i / %i\n', seed, k, N_sites);
            save(sprintf('SVDerror-seed=%i.mat', seed), 'Gs','Ls','Op','k');
            rethrow(ME)
        end
    end
end

function Di = effinv(D) % make effective inverse of diagonal matrix
    Di = diag(D);
    inds = find(Di>1e-6);
    if length(inds) ~= length(Di)
        fprintf('zero schmidt values?\n');
    end
    Di(inds) = 1./Di(inds);
    Di = full(diag(Di));
end

function tens_out = renorm(tensor)
    temp = norm(tensor(:));
%     if abs(temp-1) > 1e-6
%         fprintf('coeff tensor not normalized? norm-1 = %0.2e\n', temp-1);
%     end
    tens_out = tensor/temp;
end