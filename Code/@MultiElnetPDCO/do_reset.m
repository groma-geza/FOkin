function do_reset(self)
% sets internal values for cold restart
n = self.n;
m = self.m;
p = self.p;
n_pdco = 2*n*p;

self.x = ones(n_pdco,1)/n_pdco;
self.y = zeros(m*p,1);             
self.z = self.x;
end
