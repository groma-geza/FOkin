function optimize(self, warm)
arguments
    self
    warm(1,1) logical = false;
end
assert(self.ready, 'The FOkin object is not ready for optimization.')
self.results_ = []; % delete earliar results

p = self.p;
group_lasso = self.options_.group_lasso;
if group_lasso && (p > 1) && ~(self.optimizer.supports_group)
    error('The optimizer does not support group fit.')
end

self.optimizer.name = self.name_;

if group_lasso || (p == 1)
    self.optimizer.A = self.A;
    self.optimizer.b = self.bfit;
    self.optimizer.optimize(warm);
    x_norm = self.optimizer.result;
    self.objval = self.optimizer.objective();
    self.numiter = self.optimizer.numiter();
    self.runtime = self.optimizer.runtime();
else
    x_norm = zeros(self.n, p);
    objval = zeros(1, p);
    numiter = zeros(1, p);
    runtime = zeros(1 ,p);
    for ii = 1:p
        if self.multi_A
            self.optimizer.A = self.A{ii};
            self.optimizer.b = self.bfit{ii};
        else
            self.optimizer.A = self.A;
            self.optimizer.b = self.bfit(:,ii);
        end
        self.optimizer.optimize(warm);
        x_norm(:,ii) = self.optimizer.result;
        objval(ii) = self.optimizer.objective();
        numiter(ii) = self.optimizer.numiter();
        runtime(ii) = self.optimizer.runtime();
    end
    self.objval = objval;
    self.numiter = numiter;
    self.runtime = runtime;
end

self.x_norm = x_norm;

self.build_results('optimize')
end