function value = MPSexpectation(MPSin, Ops, k, onelocal)
%MPSexpectation calculations expectation value of operators in MPS
%   Usage:
%   value = MPSexpectation(MPSin, Ops)
%       where Ops = a cell containing N 1-local operators
%                   (where empty cell = identity)
%
%   value = MPSexpectation(MPSin, Op, site_ind, onelocal)
%       where Op = a 1- or 2-local operator
%             site_ind = index of the first site Op is acting on
%             onelocal = 1 if 1-local, 0 if 2-local (optional, default = 1)

Gs = MPSin.G;
Ls = MPSin.L;

multilocal = 0;
if iscell(Ops)
    multilocal = 1;
    if length(Ops) ~= length(Gs)
        error('number of operators do not match number of sites');
    end
end

global d

N_sites = length(Gs);

efficientTrace = @(A, B) sum(sum(conj(A).*B));


if multilocal
    G_temp = Gs;
    for k = 1:length(Ops)
        if ~isempty(Ops{k})
            G_temp = MPSOneSiteOp(G_temp, Ops{k}, k);
        end
    end
    value = MPSoverlap(Gs, Ls, G_temp, Ls);
elseif onelocal
    if k == 1
        DR = length(Ls{k});
        tensor = zeros(d, DR);
        for phys_ind = 1:d
            tensor(phys_ind, :) = Gs{1}(:,phys_ind).' * Ls{1};
        end
        value = efficientTrace(tensor, Ops*tensor);
    elseif k == N_sites
        DL = length(Ls{k-1});
        tensor = zeros(d, DL);
        for phys_ind = 1:d
            tensor(phys_ind, :) = Ls{N_sites-1} * Gs{N_sites}(:,phys_ind);
        end
        value = efficientTrace(tensor, Ops*tensor);
    else
        DL = length(Ls{k-1});
        DR = length(Ls{k});
        tensor = zeros(d, DL, DR);
        for phys_ind = 1:d
            tensor(phys_ind, :, :) = Ls{k-1} * Gs{k}(:, :, phys_ind) * Ls{k};
        end
        tensor = reshape(tensor, [d, DL*DR]);
        value = efficientTrace(tensor, Ops*tensor);
    end
else % twolocal
    if k == 1
        DM = length(Ls{1});
        DR = length(Ls{2});
%         tensor = zeros(d, d, DR);
%         for phys_ind1 = 1:d
%             for phys_ind2 = 1:d
%                 tensor(phys_ind1, phys_ind2, :) = ...
%                     Gs{1}(:, phys_ind1).' * Ls{1} * Gs{2}(:, :, phys_ind2) * Ls{2};
%                     %   1xDM                    DMxDM            DMxDR                      DRxDR
%             end
%         end
%         tensor = reshape(tensor, [d^2, DR]);

        G1 = Gs{1}.'*Ls{1};
        G2 = reshape(permute(Gs{2}, [1,3,2]), [DM*d, DR])*Ls{2};
        G2 = reshape(full(G2), [DM, d*DR]);
        tensor = reshape(G1*G2, [d^2, DR]);
        
        value = efficientTrace(tensor, Ops*tensor);
    elseif k == N_sites - 1
        DL = length(Ls{k-1});
        DM = length(Ls{k});
%         tensor = zeros(d, d, DL);
%         for phys_ind1 = 1:d
%             for phys_ind2 = 1:d
%     %             size(Ls{k-1})
%                 tensor(phys_ind1, phys_ind2, :) = ...
%                     Ls{site_ind-1} * Gs{site_ind}(:,:, phys_ind1) * Ls{site_ind} * Gs{site_ind+1}(:, phys_ind2);
%                     %  DLxDL            DLxDM                       DMxDM            DMx1
%             end
%         end
%         tensor = reshape(tensor, [d^2, DL]);
        
        G1 = Ls{k-1}*reshape(Gs{k}, [DL, DM*d]);
        G1 = reshape(permute(reshape(full(G1),[DL, DM, d]), [1,3,2]), [DL*d, DM]);
        G2 = Ls{k}*Gs{k+1};
        tensor = reshape(G1*G2, [DL, d^2]).';
        
        value = efficientTrace(tensor, Ops*tensor);
    else
        DL = length(Ls{k-1});
        DM = length(Ls{k});
        DR = length(Ls{k+1});
%         tensor = zeros(d, d, DL, DR);
%         for phys_ind1 = 1:d
%             for phys_ind2 = 1:d
%                 tensor(phys_ind1, phys_ind2, :, :) = ...
%                     Ls{site_ind-1} * Gs{site_ind}(:, :, phys_ind1) * Ls{site_ind} * Gs{site_ind+1}(:, :, phys_ind2) * Ls{site_ind+1};
%                    %   DLxDL            DLxDM          DMxDM           DMxDR               DRxDR
%             end
%         end
%         tensor = reshape(tensor, [d^2, DL*DR]);

        G1 = reshape(permute(Gs{k}, [1,3,2]), [DL, d*DM]);
        G1 = reshape(full(Ls{k-1}*G1), [DL*d, DM]);
        G2 = reshape(permute(Gs{k+1}, [1,3,2]), [DM*d, DR]);
        G2 = reshape(full(G2*Ls{k+1}), [DM, d*DR]);
        tensor = reshape(permute(reshape(G1*Ls{k}*G2, [DL, d, d, DR]), [2, 3, 1, 4]), [d^2, DL*DR]);
        value = efficientTrace(tensor, Ops*tensor);
    end
end

end