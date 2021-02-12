function QMPS = getQbprime_as_MPS(p, betas)
%getQbprime_as_MPS returns Q_{b'} as a rank-4 MPS
%
%   QMPS = getQbprime_as_MPS(p, betas)

d = 4;

%% helper

%  b_{+j}  b_{-j}
bs = [1     1;...
      1    -1;...
     -1     1;...
     -1    -1];
 
x_pL = bs(:, 1); % x_{j}
x_pR = x_pL;     % x_{j+1}
x_mL = bs(:, 2); % x_{-j}
x_mR = x_mL;     % x_{-j-1}

%%

Ls = cell(p-1,1);
Gs = cell(p,1);

myD = 4;

for r = 1:p-1
    Ls{r} = speye(myD);
end

for r = 2:p-1
    Gs{r} = zeros(myD,myD,d);

    for ind = 1:d
        b_pj = bs(ind, 1);
        b_mj = bs(ind, 2);       
        tempR = cos(betas(r)).^(1 + (b_pj * x_pR + b_mj*x_mR)/2) ...
                    .*  sin(betas(r)).^(1 - (b_pj*x_pR + b_mj*x_mR)/2) ...
                    .* (1i).^((b_mj*x_mR - b_pj*x_pR)/2);
        tempL = 0 + (x_pL == b_pj).* (x_mL == b_mj);
        
        Gs{r}(:,:,ind) = tempL * tempR.';
    end
end

Gs{1} = zeros(myD, d);
Gs{p} = zeros(myD, d);

for ind = 1:d
    b_pj = bs(ind, 1);
    b_mj = bs(ind, 2);
    Gs{1}(:, ind) = cos(betas(1)).^(1 + (b_pj*x_pR + b_mj*x_mR)/2) ...
                    .*  sin(betas(1)).^(1 - (b_pj*x_pR + b_mj*x_mR)/2) ...
                    .* (1i).^((b_mj*x_mR - b_pj*x_pR)/2);
                
    temp = cos(betas(p))^(1+ (b_pj + b_mj)/2)  ...
         * sin(betas(p))^(1- (b_pj + b_mj)/2) ...
         * (1i).^((b_mj - b_pj)/2);
    Gs{p}(:, ind) = temp * (x_pL == b_pj).* (x_mL == b_mj);
end

QMPS.G = Gs;
QMPS.L = Ls;


end