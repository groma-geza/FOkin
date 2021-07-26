function do_reset(self)
% sets internal values for cold restart
n = self.n;
p = self.p;
self.x = zeros(n, p);
self.z = zeros(n, p);
self.u = zeros(n, p);
self.rho = [];
end
