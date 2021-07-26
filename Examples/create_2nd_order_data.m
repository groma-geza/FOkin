function [time, data, txt1] = create_2nd_order_data(C, noise_sigma,...
    with_exp, show)

dt = 1.E-2;
logt=(-2:dt:4)';
time=10.^logt;
tau = time;

A = FOkin.design_matrix(time, 0, tau, 0);

k1 = 1.E-1;
k2= 1.E1;
k3 = 1.E-3;

txt1 = '2nd order';
data = C./(1 + C*k1*time);
if with_exp
    txt1 = '2nd order + exp';    
    data = data + (exp(-k2*time) + exp(-k3*time));
end

x_true = lsqminnorm(A, data);
data_rec = A*x_true;
residual = data - data_rec;
MSE = sum(residual.^2)/length(residual);
txt2 = [txt1, '  MSE = ', num2str(MSE)];

if show
    figure('Position', [50, 620, 560, 420]);
    subplot(2, 1, 1)
    semilogx(time, data, 'b', time, data_rec, 'r')
    xlabel('Time (a. u.)')
    ylabel('Signal (a. u.)')
    title(txt2)
    subplot(2, 1, 2)
    semilogx(time, residual, 'r')
    xlabel('Time (a. u.)')
    ylabel('Residual (a. u.)')
    
    figure('Position', [620, 620, 560, 420]);
    semilogx(tau, x_true, 'r.-')
    xlabel('\tau (a. u.)')
    title(txt2)
end

%% Add Gaussian noise to the data
% Set the random number generator as it was in the presented calculation
% to ensure generating identical noise.
rng(2)
randn(1, 51);
noise = noise_sigma*max(max(abs(data)))*randn(size(data));
data = data + noise;
