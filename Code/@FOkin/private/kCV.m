function [meanMSPE, user_data] = kCV(self, A, bfit, weightfit, len,...
    simple, warm)
% Calculates Mean Squared Prediction Error (MPSE) according to the
% method of k-fold cross-validation.

options = self.options_;
group_lasso = options.group_lasso;
nfold = options.cv_nfold;
num_par_workers = options.num_par_workers;
n = self.n;
p = self.p;
optimizer = self.optimizer;
MSPE = NaN(1, nfold);
accumweight = NaN(1, nfold);
x_train = zeros(n, p, nfold);

cvp = cell(1,p);
for jj = 1:p
    cvp{jj} = cvpartition(len(jj), 'kFold', nfold);
end

% Do train.
parfor (ii = 1:nfold, num_par_workers)
    % Trick to bring array and cell broadcast variables in parfor.
    dummy = {optimizer, cvp, A, bfit}; %#ok<NASGU>

    if (num_par_workers > 0)
        % Parallel pool is started.
        opt_copy = optimizer.copy(); % deep copy
    else
        % The loop runs in sequential mode.
        opt_copy = optimizer; % handle copy
    end
    
    A_train = cell(1, p);
    bfit_train = cell(1, p);
    for jj = 1:p
        train = cvp{jj}.training(ii);
        A_jj = A{jj};
        bfitjj = bfit{jj};
        A_train{jj} = A_jj(train,:);
        bfit_train{jj} = bfitjj(train);
    end
    if (p == 1) % for optimizers able to handle only this case
        opt_copy.A = A_train{1};
        opt_copy.b = bfit_train{1};
    elseif group_lasso
        opt_copy.A = A_train;
        opt_copy.b = bfit_train;
    end
    if simple
        opt_copy.optimize(warm);
        x_train(:,:,ii) = opt_copy.result;
    else
        % Warm restart is not adequate.
        x_temp = zeros(n, p);
        for jj = 1:p
            opt_copy.A = A_train{jj};
            opt_copy.b = bfit_train{jj};
            opt_copy.optimize(false);
            x_temp(:,jj) = opt_copy.result;
        end
        x_train(:,:,ii) = x_temp;
    end
end

% Do test
for ii = 1:nfold
    test_len = NaN(1, p);
    for jj = 1:p
        test_len(jj) = cvp{jj}.TestSize(ii);
    end
    residual_test = NaN(sum(test_len), 1);
    weightfit_test = NaN(sum(test_len), 1);
    
    ptr = 1;
    for jj = 1:p
        test_jj = cvp{jj}.test(ii);
        A_jj = A{jj};
        fit_test_jj = A_jj(test_jj,:)*x_train(:,jj,ii);
        bfit_test_jj = bfit{jj}(test_jj);
        weightfit_test_jj = weightfit{jj}(test_jj);
        
        % Compensate the effect of the method applied for acheaving a
        % weighted fit in self.rebuild()
        residual_test_jj = (bfit_test_jj - fit_test_jj)./...
            sqrt(weightfit_test_jj);
        
        % The order of the residuals is left random, we need only
        % the weighted average of them.
        new_ptr = ptr + test_len(jj);
        residual_test(ptr:(new_ptr-1)) = residual_test_jj;
        weightfit_test(ptr:(new_ptr-1)) = weightfit_test_jj;
        ptr = new_ptr;       
    end
    
    % Mean Squared Prediction Error
    MSPE(ii) = wmse(residual_test, weightfit_test);
    accumweight(ii) = sum(weightfit_test);
end

% The cross validation is characterized by the weighted average and
% stardard deviation of mean square errors corresponding to the individual
% tests.

meanMSPE = sum(MSPE.*accumweight)/sum(accumweight);
stdMSPE = sqrt(sum(((MSPE - meanMSPE).^2).*accumweight)/sum(accumweight)...
    /(nfold - 1));

x_train_avr = sum(x_train, 3)/nfold;
% Make x_train_avr available for self.results().
self.x_norm = x_train_avr;

user_data.stdMSPE = stdMSPE;
end


