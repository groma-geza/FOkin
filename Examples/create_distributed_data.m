function [time, data] = create_distributed_data(noise_sigma, show)

%% Set Arrhenius equation with distributed activation energy
E0=200;
s=35;

E=(0:400)';
g=exp(-(E-E0).^2/2/s/s)/s/sqrt(2*pi);
k=exp(-E/50);

%% Calculate the true distribution of tau
g_ext = [0; g; 0];
tau = 1./k;
tau_ext = [1.E-1; tau; 1.E4];

%% Calculate a signal proportional to the concentration 
logt=(-1:.1:4)';
time=10.^logt;

A = FOkin.design_matrix(time, 0, tau, 0);
b_dist = A*g;
b = exp(-exp(-E0/50)*time);

data = b_dist;

%% Add Gaussian noise to the data
% Set the random number generator as it was in the presented calculation
% to ensure generating identical noise.
rng(2)
randn(1, 51);
noise = noise_sigma*max(max(abs(data)))*randn(size(data));
data = data + noise;
rng('shuffle') % shuffle the random number generator
               % for further calculations
               
if show
    figure % S3 Fig A
    yyaxis left
    plot(E, g); grid on
    xlabel('Activation Energy (a. u.)');
    ylabel('Distribution')
    hold on
    yyaxis right
    plot(E,k);
    grid on
    ylabel('Rate constant (a. u.)')

    figure % S3 Fig B, Fig 11 red line
    semilogx(tau_ext,g_ext);
    grid on
    xlabel('\tau (a. u.)');
    ylabel('Amplitude')

    figure % Fig 8
    semilogx(time, b_dist, time, b)
    grid on
    xlabel('Time (a. u.)');
    ylabel('Signal')
    legend('distributed','discrete')
end
