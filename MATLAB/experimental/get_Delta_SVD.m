function [myU, myV, myS] = get_Delta_SVD(aList, gammas)
%get_Delta_SVD returns the (unnormalized) Singular-Value Decomposition of
%   Delta_{aPrime, bPrime} matrix
%
%   Example:
%       
%   Let
%       aPrime = unstar_fun(aList);
%       Dmat = Delta_fun(aPrime, aPrime, gammas);
%       [myU, myV, myS] = get_Delta_SVD(aList, gammas);
%
%   Then (up to floating point errors):
%       norm(myU* diag(myS) * myV' - Dmat, 'fro') == 0

    p = size(aList,2) / 2;
    
    fp = @(ind) (aList(:, ind) + aList(:, 2*p+1-ind))/2;
    fm = @(ind) (aList(:, ind) - aList(:, 2*p+1-ind))/2;

    ijs = nchoosek(1:p,2);

    myU = zeros(4^p, p*(p-1));
    myV = myU;
    myS = zeros(p*(p-1),1);


    for ij_ind = 1:size(ijs,1)
        ind = ijs(ij_ind,1);
        ind2 = ijs(ij_ind,2);

    %     fprintf('(%d, %d), (%d, %d) = %d\n', ind, 2*p+1-ind, ind2, 2*p+1-ind2, ij_ind);

        myU(:, ij_ind) = fm(ind) .* fp(ind2);
        myV(:, ij_ind) = fp(ind) .* fm(ind2);
        for ind3 = ind2+1:p
           myU(:, ij_ind) = myU(:, ij_ind) .* fp(ind3).^2;
        end
        myS(ij_ind) = -8*gammas(ind)*gammas(ind2);

        if ind2-ind >= 2

            myind = p*(p-1)/2 + ij_ind;
            myU(:, myind) = fp(ind) .* fp(ind2);

            for ind3 = ind2+1:p
               myU(:, myind) = myU(:, myind) .* fp(ind3).^2;
            end

            temp = ones(4^p,1);
            for ind3 = ind+1:ind2-1
                temp = temp .* fp(ind3).^2;
            end
            myU(:, myind) = myU(:, myind).* (1-temp);

            myV(:, myind) = fm(ind) .* fm(ind2);

            myS(myind) = -8*gammas(ind)*gammas(ind2);
        end
    end


    Ivalid = myS ~= 0;
    assert(nnz(Ivalid) == (p-1)^2);
    myU = myU(:, Ivalid);
    myV = myV(:, Ivalid);
    myS = myS(Ivalid);

end