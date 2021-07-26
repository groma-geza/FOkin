function [meanMSPE, user_data] = RCVnv(self, A, bfit, weightfit, len,...
    simple, warm)
% Calculates Mean Squared Prediction Error (MPSE) according to the
% restricted leave-n_v-out cross-validation method desribed as Algoritm 2
% in Y. Feng and % Y. Yu "The restricted consistency property of 
% leave-nv-out cross-validation for high-dimensional variable selection" 
% Statistica Sinica, 29, 1607-1630 (2019), referred below as F&YA2

options = self.options_;
group_lasso = options.group_lasso;
nrep = options.cv_nrep;
tol = options.lsqminnorm_tol;
nc_rel = options.nc_rel;
num_par_workers = options.num_par_workers;
n = self.n;
p = self.p;
optimizer = self.optimizer;
MSPE = NaN(1, nrep);
accumweight = NaN(1, nrep);

test_len = NaN(1, p);
for jj = 1:p
    test_len(jj) = len(jj) - round(nc_rel*len(jj));
end

A_supp = cell(1, p);

% F&YA2 S1
%Instead of building a solution path, the lambda and alpha
% penalization parameters are choosen by the 'bayesopt' function.
% Obtain the support of the solution by solving the penalized problem.
if (p == 1) % for optimizers able to handle only this case
    optimizer.A = A{1};
    optimizer.b = bfit{1};
elseif group_lasso
    optimizer.A = A;
    optimizer.b = bfit;
end
if simple
    optimizer.optimize(warm);
    x_whole = optimizer.result;
else
    % Warm restart is not adequate.
    x_whole = NaN(n, p);
    for jj = 1:p
        optimizer.A = A{jj};
        optimizer.b = bfit{jj};
        optimizer.optimize(false);
        x_whole(:,jj) = optimizer.result;
    end
end
support = (abs(x_whole) > options.nonzero_limit); % \alpha in F&YA2

ns = sum(support, 1);

if group_lasso && any(diff(sum(support)))
    % Something went wrong in the optimization.
    error('The obtained supports are not identical for all groups.')
end

for jj = 1:p
    A_supp{jj} = A{jj}(:, support(:,jj));  
end

x_train_supp_avr = zeros(1, sum(ns));

% Do repeated Maximum Likelihood train and test only on the support. 

% for ii = 1:nrep
 parfor (ii = 1:nrep, num_par_workers)
    % Trick to bring array and cell broadcast variables in parfor.
    dummy = {len, A_supp, bfit, weightfit}; %#ok<NASGU>
    
    residual_test = NaN(sum(test_len), 1);
    weightfit_test = NaN(sum(test_len), 1);
    x_train = NaN(1, sum(ns));
    ptr1 = 1;
    ptr2 = 1;
    for jj = 1:p
        % F&YA2 S2
        % Create random split separately in each group
        n_jj = len(jj);
        nc_jj = round(nc_rel*n_jj);
        idx = randperm(n_jj);
        train_jj = idx(1:nc_jj);
        test_jj = idx((nc_jj+1):end);        
        A_supp_jj = A_supp{jj};
        
        if ~isempty(A_supp_jj)
            
            % F&YA2 S3 (a)
            % Do train
            A_train_jj = A_supp_jj(train_jj,:);
            bfit_train_jj = bfit{jj}(train_jj);
            
            % For the Gaussian likelihood the solution of Eq. 3.1 in F&YA2
            % corresponds to the solution of a linear least square problem
            % on the support:
            % bfit_train_jj = A_train_jj*x_train_jj
            %
            % Apply lsqminnorm which minimizes the the norm of the residual
            % if the problem is underdetermined. Optionally specify a
            % tolerance for the underlining complete orthogonal 
            % decomposition algorithm.

            if isempty(tol)
                x_train_jj = lsqminnorm(A_train_jj, bfit_train_jj);
            else
                x_train_jj = lsqminnorm(A_train_jj, bfit_train_jj, tol);
            end
           
        else
            [~, col] = size(A_supp_jj);
            x_train_jj = zeros(col, 1);
        end
        
        % F&YA2 S3 (b)
        % Do test
        A_test_jj = A_supp_jj(test_jj,:);
        bfit_test_jj = bfit{jj}(test_jj);
        fit_test_jj = A_test_jj*x_train_jj;
        weightfit_test_jj = weightfit{jj}(test_jj);
        
        % Compensate the effect of the method applied for acheaving a
        % weighted fit in self.rebuild()
        residual_test_jj = (bfit_test_jj - fit_test_jj)./...
            sqrt(weightfit_test_jj);

        new_ptr = ptr1 + test_len(jj);
        residual_test(ptr1:(new_ptr-1)) = residual_test_jj;
        weightfit_test(ptr1:(new_ptr-1)) = weightfit_test_jj;
        ptr1 = new_ptr;
        
        new_ptr = ptr2 + ns(jj);
        x_train(ptr2:(new_ptr-1)) = x_train_jj;
        ptr2 = new_ptr;
    end

    % Mean Squared Prediction Error
    MSPE(ii) = wmse(residual_test, weightfit_test);
    accumweight(ii) = sum(weightfit_test);
    
   x_train_supp_avr = x_train_supp_avr + x_train/nrep; 
end

% F&YA2 S4
meanMSPE = sum(MSPE.*accumweight)/sum(accumweight);
stdMSPE = sqrt(sum(((MSPE - meanMSPE).^2).*accumweight)/sum(accumweight)...
    /(nrep - 1));

x_train_avr = zeros(size(x_whole));
ptr = 1;
for jj = 1:p
    new_ptr = ptr + ns(jj);
    x_train_avr(support(:,jj),jj) = (x_train_supp_avr(ptr:(new_ptr-1)))';
    ptr = new_ptr;
end

% Make x_train_avr available for self.results().
self.x_norm = x_train_avr;

user_data.stdMSPE = stdMSPE;
user_data.support_size = ns;
end

