function Delta = Delta_fun(aConfigs, bConfigs, gammas)
%Delta_fun a "vectorized" function to calculate \Delta_{ab} for a list of a
%       and b bit-configurations, and specified gammas
%
%   \Delta_{a,b} = (Phi^2_{a star * b} - Phi^2_{a * b} ) / 2
%
%   Usage:
%     Delta = Delta_fun(aConfigs, bConfigs, gammas)
%
%   Suppose aConfigs is an La x 2p matrix consisting of La rows of configurations
%           bConfigs is an Lb x 2p matrix consisting of Lb rows of configurations
%   Then Delta is an La x Lb matrix whose elements are Delta_{ab} with
%       a, b from corresponding row of aConfigs and bConfigs
%
    Delta = zeros(size(aConfigs,1), size(bConfigs,1));
    
    if size(aConfigs, 1) > size(bConfigs, 1)
        for ind = 1:size(bConfigs, 1)
            Delta(:, ind) = ( phi2_fun(star_fun(abar_fun(aConfigs) .* bConfigs(ind,:)), gammas) - ...
                              phi2_fun(star_fun(aConfigs .* bConfigs(ind,:)), gammas) ) / 2;
        end
    else
        for ind = 1:size(aConfigs, 1)
            Delta(ind, :) = ( phi2_fun(star_fun(abar_fun(aConfigs(ind,:)) .* bConfigs), gammas) - ...
                              phi2_fun(star_fun(aConfigs(ind,:) .* bConfigs), gammas) ) / 2;
        end
    end
end