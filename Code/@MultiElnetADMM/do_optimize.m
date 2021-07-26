function do_optimize(self)

% Solves the group-lasso problems by Alternating Direction Method of
% Multipliers described in the paper
% Boyd, S.; Parikh, N.; Chu, E.; Peleato, B.; Eckstein, J.,
% Distributed Optimization and Statistical Learning via the Alternating 
% Direction Method of Multipliers. Foundations and Trends® in Machine 
% Learning 2011, 3, 1-122.

% This code is a modified version of the group_lasso.m function available in
% the collection 'MATLAB scripts for alternating direction method of
% multipliers', https://web.stanford.edu/~boyd/papers/admm/, written by the
% above authors.

tstart = tic;
% local variables
% The iteration runs faster with local variables than with class variables!
opt = self.options;
rho = opt.RHO0;
alpha_relax = opt.ALPHA_RELAX;
ABSTOL   = opt.ABSTOL;
RELTOL   = opt.RELTOL;
MAX_ITER = opt.MAX_ITER;
DIAGNOSE  = opt.DIAGNOSE;

self.h_diag_fig = [];
m = self.m;
n = self.n;
p = self.p;
lambda1 = self.lambda1;
lambda2 = self.lambda2;
multi_A = self.multi_A;
A = self.A;
b = self.b;
Atb = self.Atb;

% Check if warm restart is possible
assert((size(self.x, 1) == n), ['Optimization with warm restart is ',...
    'not possible due to incompatible sizes.']);

x = self.x;
u = self.u;
if ~isempty(self.rho) % warm restart
    u = u*self.rho/rho;
end
z = self.z;

%% ADMM solver

factorize = true;
numiter = 0;
history.DIAGNOSE = DIAGNOSE;

l = 0;
longcount = 0;

if multi_A
    L_multi = cell(1, p);
    U_multi = cell(1, p);
end

for k = 1:MAX_ITER
    numiter = numiter + 1;
    
    % x-update
    rho_lambda = rho + lambda2;
    if multi_A
        for ii = 1:p
            Aii = A{ii};
            if factorize
                [L_multi{ii}, U_multi{ii}] = factor(Aii, rho_lambda);
            end
            L = L_multi{ii};
            U = U_multi{ii};
            q = Atb{ii} + rho*(z(:,ii) - u(:,ii)); % temporary
            if( m >= n )    % if skinny
                x(:,ii) = update_skinny(q, U, L);
            else            % if fat
                x(:,ii) = update_fat(q, rho_lambda, Aii, U, L);                
            end
        end
    else
        if factorize
            [L, U] = factor(A, rho_lambda);
        end
        q = Atb + rho*(z - u); % temporary
        if( m >= n )    % if skinny
            x = update_skinny(q, U, L);
        else            % if fat
            x = update_fat(q, rho_lambda, A, U, L);
        end  
    end
    factorize = false;
    
    % z-update
    zold = z;
    x_hat_matr = alpha_relax*x + (1-alpha_relax)*zold;
    z = shrinkage(x_hat_matr + u, lambda1/rho);
    
    % u-update
    u = u + (x_hat_matr - z);
    
    % diagnostics
    r_norm  = norm(x(:) - z(:));
    s_norm  = norm(-rho*(z(:) - zold(:)));
    eps_pri = sqrt(n*p)*ABSTOL + RELTOL*max(norm(x(:)), norm(-z(:)));
    eps_dual= sqrt(n*p)*ABSTOL + RELTOL*norm(rho*u(:));
    history.r_norm(k) = r_norm;
    history.s_norm(k) = s_norm;
    history.eps_pri(k) = eps_pri;
    history.eps_dual(k) = eps_dual;
    
    % reporting
    if DIAGNOSE > 0
        history.objval(k) = self.objective(z);
        history.objval_ADMM(k) = objective_ADMM(multi_A, A, b,...
            lambda1, lambda2, x, z);
        
        history.nonzero(k) = numel(find(sum(abs(z), 2) > 0));
        
        
    end
    
%     % rho-update
    step = 200;
    longsection = 10;

    if (mod(numiter, step) == 0)
        l = l + 1;
        longcount = longcount + 1;
        factorize = false;
        ratio = 1;
        s_norm_mean = mean(history.s_norm((k-step+1):k));
        if all(history.r_norm((k-step+1):k) > 0) && all(history.s_norm((k-step+1):k) > 0)
            eps_dual_mean = mean(history.eps_dual((k-step+1):k));
            dist_dual = log10(s_norm_mean/eps_dual_mean);
            if dist_dual < 0  
                ratio = 10;
                longcount = 1;
            end
            if longcount == longsection
                if rho > 1.E-6
                    ratio = 1/3; % avoid oscillation with e.g ratio = 1/10;
                else
                    ratio = 10;
                end
                longcount = 1;
            end
        elseif s_norm_mean == 0
            ratio = 10;
            longcount = 1;
        end
        
        if ratio ~= 1
            factorize = true;
            rho = rho*ratio;
            u = u/ratio;
        end
        
        if DIAGNOSE > 0
            history.rho_numiter(l) = numiter;
            history.rho_ratio(l) = ratio;
            history.rho_value(l) = rho;
            self.numiter_ = numiter;
            self.history = history;
        end
        if DIAGNOSE > 1
            self.diagnose()
        end
    end
    
    % termination checks
    if ((r_norm < eps_pri) && (s_norm < eps_dual))
        self.numiter_ = numiter;
        self.history = history;
        if DIAGNOSE > 1
            self.diagnose()
        end
        break;
    end
end

if numiter == MAX_ITER
    error('Maximum number of iterations is reached.')
end

% put the state of the iteration back to the class variables
self.runtime_ = toc(tstart);
self.x = x;
self.u = u;
self.z = z;
self.rho = rho;
self.numiter_ = numiter;
self.history = history;
end

%% local functions

function objval = objective_ADMM(multi_A, A, b, lambda1, lambda2, x, z)
% objective function on x and z, characterizing the ADMM iterations
% During the iteration x is updated to minimize norm1
% z is updated later, hence its value is better for norm2 and norm3
% At the end of iteration x and z are almost equal

    if multi_A
        p = length(b);
        q = zeros(1, p);
        for ii = 1:p
            q(ii) = sum((A{ii}*x(:,ii) - b{ii}).^2);
        end
    else
        q = sum((A*x - b).^2, 1);
    end
    norm1 = sum(q); % square of the L2 norm of the residual
    norm2 = sum(sqrt(sum(z.^2, 2))); % group-lasso penalty
    norm3 = sum(sum(z.^2)); % L2 penalty
    objval = 1/2*norm1 + lambda1*norm2 + 1/2*lambda2*norm3;
end

%%
function z = shrinkage(x, kappa)
    norm = sqrt(sum(x.^2, 2)); % L2 norm of the rows of x
    z = max((1 - kappa./norm), 0).*x;
end

%%
function [L, U] = factor(A, rho_lambda)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho_lambda*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho_lambda*(A*A'), 'lower' );
    end
    
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

%%
function x = update_skinny(q, U, L)
    x= U \ (L \ q);
end

%%
function x = update_fat(q, rho_lambda, A, U, L)
    x = q/rho_lambda - (A'*(U \ ( L \ (A*q) )))/rho_lambda^2;
end
