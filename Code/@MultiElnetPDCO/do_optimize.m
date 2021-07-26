function do_optimize(self)

n = self.n;
m = self.m;
n_pdco = 2*n;

pdObj = self.lambda1*ones(n_pdco,1); % phi(x) is linear (L1 penalty)
pdMat = [self.A, -self.A]; % positive and negative part
bl = 1.E-8*ones(n_pdco,1);
bu = Inf*ones(n_pdco,1);
mind1 = 0; % Minimum value of d1
d1 = max(mind1, sqrt(self.lambda2)); % L2 penalty
d2 = 1;
xsize = 1;                 
zsize = 1;                 

% Check if warm restart is possible
if length(self.x) ~= n_pdco
    error(['Optimization with warm restart is not possible',...
        ' due to incompatible sizes']);
end
if length(self.y) ~= m
    self.y = zeros(m,1); %partially warm restart
end

tstart = tic;
[self.x, self.y, self.z, inform, PDitns] = ...
    pdco(pdObj, pdMat, self.b, bl, bu, d1, d2, self.options, ...
    self.x, self.y, self.z, xsize, zsize);
self.runtime_ = toc(tstart);
self.numiter_ = PDitns;

switch inform
    case 0
        % normal exit
    case 1
        error('In PDCO too many iterations were required.')
    case 2
        error('In PDCO the linesearch failed too often.')
    case 3
        error('In PDCO the step lengths became too small.')
    case 4
        error('In PDCO Cholesky said ADDA was not positive definite.')
end
end

