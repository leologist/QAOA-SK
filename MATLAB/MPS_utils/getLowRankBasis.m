function [myUs, myVs, myLs] = getLowRankBasis(p)
%getLowRankBasis returns the low_rank basis vectors as rank-1 or 2 MPS
%
%   [myUs, myVs, myLs] = getLowRankBasis(p)

d = 4;

%% helper
b_p = [1,1,-1,-1];
b_m = [1,-1,1,-1];

L0 = cell(p-1,1);
for jj = 1:p-1
    L0{jj} = 1;
end


%% get Us and Vs

myUs((p-1)^2).G = [];
myUs((p-1)^2).L = [];
myUs((p-1)^2).label = [];
myUs((p-1)^2).type = 'U';
myVs = myUs;

jk_pairs = nchoosek(1:p, 2);
for ii = 1:size(jk_pairs,1)
    j = jk_pairs(ii, 1);
    k = jk_pairs(ii, 2);

    [G, L] = U_as_MPS(j, k);
    myUs(ii).G = G;
    myUs(ii).L = L;
    myUs(ii).label = [j, k, +1];
    myUs(ii).type = 'U';

    [G, L] = V_as_MPS(j, k, +1);
    myVs(ii).G = G;
    myVs(ii).L = L;
    myVs(ii).label = [j, k, +1];
    myVs(ii).type = 'V';
end

base_ind = size(jk_pairs,1);

Idelete = diff(jk_pairs,[],2) == 1;
jk_pairs(Idelete,:) = [];

for ii = 1:size(jk_pairs,1)
    j = jk_pairs(ii, 1);
    k = jk_pairs(ii, 2);

    [G, L] = Utilde_as_MPS(j, k);
    myUs(base_ind + ii).G = G;
    myUs(base_ind + ii).L = L;
    myUs(base_ind + ii).label = [j, k, -1];
    myUs(base_ind + ii).type = 'U';

    [G, L] = V_as_MPS(j, k, -1);
    myVs(base_ind + ii).G = G;
    myVs(base_ind + ii).L = L;
    myVs(base_ind + ii).label = [j, k, -1];
    myVs(base_ind + ii).type = 'V';
end

%% get Ls

myLs(2*p).G = [];
myLs(2*p).L = [];
myLs(2*p).label = [];
myLs(2*p).type = 'L_±';

for ii = 1:p
    [G, L] = L_as_MPS(ii, +1);
    myLs(ii).G = G;
    myLs(ii).L = L;
    myLs(ii).label = [ii, +1];
    myLs(ii).type = 'L_±';
    
    [G, L] = L_as_MPS(ii, -1);
    myLs(p+ii).G = G;
    myLs(p+ii).L = L;
    myLs(p+ii).label = [ii, -1];
    myLs(p+ii).type = 'L_±';
end



%% helper functions

    function [G, L] = U_as_MPS(j, k)
        if j > k
            l = j; j = k; k = l;
        end
        
        G = cell(p,1);
        for ind = 1:p
            if ind == j
                G{ind} = (b_p - b_m)/2;
            elseif ind == k
                G{ind} = (b_p + b_m)/2;
            elseif ind >= k+1
                G{ind} = (b_p + b_m).^2/4;
            else
                G{ind} = ones(1,d);
            end

            if ind >= 2 && ind <= p-1
                G{ind} = reshape(G{ind}, 1, 1, d);
            else
                G{ind} = reshape(G{ind}, 1, d);
            end
        end
        L = L0;
    end



    function [G, L] = Utilde_as_MPS(j, k)
        if j > k
            l = j; j = k; k = l;
        end
        if abs(k-j) <= 1
            error('j=%d, k=%d are invalide for Utilde', j, k);
        end
        
        L = L0; % all one scalars

        G = cell(p,1);
        for ind = 1:p
            if ind == j
                G{ind} = zeros(1,2, d);
                G{ind}(1,1,:) = (b_p + b_m)/2;
                G{ind}(1,2,:) = (b_p + b_m)/2;
                L{ind} = speye(2);
            elseif j < ind && ind < k
                G{ind} = zeros(2,2,d);
                G{ind}(1,1,:) = 1;
                if ind == j+1
                    G{ind}(2,2,:) = -(b_p + b_m).^2/4;
                else
                    G{ind}(2,2,:) = (b_p + b_m).^2/4;
                end
                L{ind} = speye(2);
            elseif ind == k
                G{ind} = zeros(2,1,d);
                G{ind}(1,1,:) = (b_p + b_m)/2;
                G{ind}(2,1,:) = (b_p + b_m)/2;
            elseif ind >= k+1
                G{ind} = zeros(1,1,d);
                G{ind}(1,1,:) = (b_p + b_m).^2/4;
            else % ind < j
                G{ind} = ones(1,1,d);
            end

            if ind == 1 || ind == p
                G{ind} = reshape(G{ind}, [], d);
            end
        end
    end

    function [G, L] = V_as_MPS(j, k, sign)
        if j > k
            l = j; j = k; k = l;
        end
        
        G = cell(p,1);
        for ind = 1:p
            if ind == j
                G{ind} = (b_p + sign * b_m)/2;
            elseif ind == k
                G{ind} = (b_p - b_m)/2;
            else
                G{ind} = ones(1,d);
            end

            if ind >= 2 && ind <= p-1
                G{ind} = reshape(G{ind}, 1, 1, d);
            else
                G{ind} = reshape(G{ind}, 1, d);
            end
        end
        L = L0;
    end

    function [G, L] = L_as_MPS(j, sign)
        L = L0;
        G = cell(p,1);
        for ind = 1:p
            if ind == j
                G{ind} = (b_p + sign * b_m)/2;
            else
                G{ind} = ones(1,d);
            end

            if ind >= 2 && ind <= p-1
                G{ind} = reshape(G{ind}, 1, 1, d);
            else
                G{ind} = reshape(G{ind}, 1, d);
            end
        end
    end

end