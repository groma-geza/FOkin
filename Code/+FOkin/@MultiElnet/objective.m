function objval = objective(self, x)
% objective function of the multi elnet problem
% input: x          - the matrix of the unknowns to be optimized. This  
%                     optional input argument can be useful e.g. during an
%                     itaration process when the class properties are not
%                     properly updated for the calculation of 
%                     self.result, normally used to obtain x.
%
% output: objval    - value of the objective function
arguments
    self
    x double = self.result;
end

p = self.p;

if isempty(x)
    objval = [];
    return
end

if self.multi_A
    q = zeros(1, p);
    for ii = 1:p
        q(ii) = sum((self.A{ii}*x(:,ii) - self.b{ii}).^2);
    end
else
    q = sum((self.A*x - self.b).^2, 1);
end
norm1 = sum(q); % square of the L2 norm of the residual
norm2 = sum(sqrt(sum(x.^2, 2))); % group-lasso penalty
norm3 = sum(sum(x.^2)); % L2 penalty
objval = 1/2*norm1 + self.lambda1*norm2 + 1/2*self.lambda2*norm3;
end