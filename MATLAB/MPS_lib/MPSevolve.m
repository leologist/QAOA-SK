function [values, jumps, err_trunc] = ...
    MPSevolve(MPS, H_2bdy, Lijs, observables, N_steps, dT, dp, samp)
%MPSEvolve: evolve MPS according to 2-body Hamiltonians and quantum jump
%           operators, in open boundary condition
%
%   Usage: [values, jumps, err_trunc] = ...
%   MPSEvolve(G_0, L_0, H_2bdy, Lijs, observables, N_steps, dT, dp, samp)
%   Output:
%       values = observables values (or state history, if observables = {})
%       jumps = number of quantum jumps over time
%       err_trunc = MPS trunction error over time
%   
%   Input:
%       G_0, L_0 = initial state given in MPS form (cell arrays of matrices)
%       H_2bdy = two-body Hamiltonians (can be a cell array of Hamiltonians
%                corresponding to pairs of sites)
%       Lijs = cell array of quantum jump operators
%       observables = cell array of observables (each element can be a cell
%                     array of operators corresponding to each site, or a
%                     function with two inputs (G,L)). If empty, observable
%                     values are replaced with state history
%       N_steps = number of time steps
%       dT = size of time step
%   
%   Optional Input:
%       dp (default 0.1) = threshold for quantum jump probability (can be [])
%       samp (default 1) = sampling period (returns observable values and
%                          other output every samp*dT)
%

G_0 = MPS.G;
L_0 = MPS.L;
%% parameters and debugging
if ~exist('dp','var') || isempty(dp)
    dp = 0.1;
end

if ~exist('samp', 'var')
    samp = 1;
end


%% housekeeping
N_sites = length(MPS.G);
num_obs = length(observables);
N_values = floor(N_steps/samp);

num_loss_ops_per_site = length(Lijs);
num_loss_ops = (N_sites - 1)*num_loss_ops_per_site;

LLs = cell(length(Lijs), 1);
for ind = 1:length(Lijs)
    LLs{ind} = Lijs{ind}'*Lijs{ind};
end

if ~iscell(H_2bdy) % translation-invariant Hamiltonian
    [eigvec, eigval] = eig(full(H_2bdy));
    eigval = diag(eigval);
    eigvec_inv = inv(eigvec);

    U_2bdy_dT2 = expm(-1i*H_2bdy*dT/2);
    U_2bdy_dT = expm(-1i*H_2bdy*dT);
    
    ED_2bdy = {U_2bdy_dT2, U_2bdy_dT, eigvec, eigval, eigvec_inv};
else
    ED_2bdy = cell(N_sites-1, 5);
    for ind = 1:N_sites-1
        ED_2bdy{ind, 1} = expm(-1i*H_2bdy{ind}*dT/2);
        ED_2bdy{ind, 2} = expm(-1i*H_2bdy{ind}*dT);        
        [eigvec, eigval] = eig(full(H_2bdy{ind}));
        eigval = diag(eigval);
        eigvec_inv = inv(eigvec);
        
        ED_2bdy(ind,3:5) = {eigvec, eigval, eigvec_inv};
    end
end

%% start evolution
jump_cache = 0;
T_final = N_steps*dT;
MPS_t = MPS;

if num_obs > 0
    values = zeros(N_values+1, num_obs);
    values(1, :) = getObsValues(MPS);
else
    values = cell(N_values+1,1);
    values{1} = {MPS};
end

jumps = zeros(N_values+1, 1);
err_trunc = zeros(N_values+1, 1);
fid_trunc = 1;

t = 0;
T_index = 0;
array_index = 1;

time_marker_increment = 60;
next_time_marker = time_marker_increment;
tic;

    while t < T_final
        dpdt_k = zeros(num_loss_ops, 1);
        for ind = 1:num_loss_ops
            site_ind = ceil(ind/num_loss_ops_per_site);
            op_ind = ind - (site_ind-1)*num_loss_ops_per_site;
            temp = MPSexpectation(MPS_t, LLs{op_ind}, site_ind, 0);
            dpdt_k(ind) = abs(temp);
        end
        dpdt = sum(dpdt_k);

        % choose dt
        tryRecordObservables = 0;
        if dpdt > 0
            dt = dp/dpdt;
            if floor((t+dt)/dT) > T_index
                dt = (T_index + 1)*dT - t;
                tryRecordObservables = 1;
            end
        else
            dt = (T_index + 1)*dT - t;
            tryRecordObservables = 1;
        end

        % evolve
        
        temp = rand; % random number between 0 and 1
        ind_k = find(temp <= (cumsum(dpdt_k)*dt), 1, 'first');

        if isempty(ind_k) % evolution by non-Hermitian Hamiltonian
            if abs(dt-dT) > 1e-14
%                 U_2bdy2 = eigvec * diag(exp(-1i*eigval*dt/2)) * eigvec_inv;
%                 U_2bdy = eigvec * diag(exp(-1i*eigval*dt)) * eigvec_inv;
                if ~iscell(H_2bdy) % translation-invariant Hamiltonian
                    U_2bdys = cell(1,2);
                    U_2bdys{1} = ED_2bdy{3}*diag(-1i*ED_2bdy{4}*dt/2) * ED_2bdy{5};
                    U_2bdys{2} = ED_2bdy{3}*diag(-1i*ED_2bdy{4}*dt) * ED_2bdy{5};
                else
                    U_2bdys = cell(N_sites-1,2);
                    for ind = 1:N_sites-1
                        U_2bdys{ind,1} = ED_2bdy{ind,3}*diag(-1i*ED_2bdy{ind,4}*dt/2) * ED_2bdy{ind,5};
                        U_2bdys{ind,2} = ED_2bdy{ind,3}*diag(-1i*ED_2bdy{ind,4}*dt) * ED_2bdy{ind,5};
                    end
                end
            else
%                 U_2bdy2 = U_2bdy_dT2;
%                 U_2bdy = U_2bdy_dT;
                U_2bdys = ED_2bdy(:,1:2);
            end

            for ind = 2:2:N_sites-1 %even
                [MPS_t, err] = MPSTwoSiteOp(MPS_t, getUnitary(ind, 1), ind, 1);
                fid_trunc = fid_trunc*(1-err);
            end
            for ind = 1:2:N_sites-1 %odd
                [MPS_t, err] = MPSTwoSiteOp(MPS_t, getUnitary(ind, 2), ind, 1);
                fid_trunc = fid_trunc*(1-err);
            end
            for ind = 2:2:N_sites-1 %even
                [MPS_t, err] = MPSTwoSiteOp(MPS_t, getUnitary(ind, 1), ind, 1);
                fid_trunc = fid_trunc*(1-err);
            end
        else % quantum jump
            site_ind = ceil(ind_k/num_loss_ops_per_site);
            op_ind = ind_k - (site_ind-1)*num_loss_ops_per_site;
            [MPS_t, err] = MPSTwoSiteOp(MPS_t, Lijs{op_ind}, site_ind, 1);
            fid_trunc = fid_trunc*(1-err);
            jump_cache = jump_cache + 1;
        end

        MPS_t = MPScanonicalize(MPS_t);

        if tryRecordObservables
            T_index = T_index + 1; % increment T_index
            
            if mod(T_index, samp) == 0
                array_index = array_index + 1; % get array_index for observables
                err_trunc(array_index) = 1-fid_trunc;
                if num_obs > 0
                    values(array_index, :) = getObsValues(MPS_t);
                else
                    values{array_index} = {MPS_t};
                end
                jumps(array_index) = jump_cache;
                jump_cache = 0;

                elapsedTime = toc();
                if elapsedTime > next_time_marker
                    fprintf('Evolution %0.1f percent complete after %0.1f min...\n', T_index/N_steps*100, elapsedTime/60);
                    next_time_marker = next_time_marker + time_marker_increment;
                end
            end

        end

        t = t + dt;
    end


    function tempvalues = getObsValues(MPS)
        tempvalues = zeros(1, num_obs);
        for obs_ind = 1:num_obs
            if isa(observables{obs_ind}, 'function_handle')
                tempvalues(obs_ind) = observables{obs_ind}(MPS);
            else                
                tempvalues(obs_ind) = MPSexpectation(MPS, observables{obs_ind});
            end
        end
    end

    function U = getUnitary(site_ind, dT_ind)
        if ~iscell(H_2bdy) % translation-invariant Hamiltonian
            U = U_2bdys{dT_ind};
        else
            U = U_2bdys{site_ind,dT_ind};
        end
    end
end