function raw = create_raw_bR_data()

%% Load spectral data
load('bR_spectral_data.mat', 'raw_wavelength', 'sp_abs', 'lamp')
sp_dif = sp_abs(:,1:end-1) - sp_abs(:,end);
% Assign the spectra to the components
sp_dif_assign = [1 1 2 3 3 4 4 4 5 5 6]; 
% (K1,K2), L1, (L2,L3), (M1,M2,M3), (N1,N2), bR

%% Define equidistant log scale time points
log_time=(-7:.115:-1.3);
time=10.^log_time;

%% Load microscopic rate constant data
% The photocycle model is defined by the logarithm of the nondiagonal
% microscopic rate constants.
log_nondiag_micro_rates = dlmread('bR_rate_data.txt');
nondiag_micro_rates = 10.^log_nondiag_micro_rates;

%% Built the complate K matrix by calculatting the diagonal rates
K = (-diag(sum(nondiag_micro_rates, 1)) + nondiag_micro_rates); % S1 Table
[n_inter, ~] = size(K) ; % number of intermediates with the
                                   % initial state
                                   
%% Solve the eigenvalue problem of K                          
[eig_vect, eig_val] = eig(K,'nobalance');
eig_val=diag(eig_val); 
true_tau = -1./eig_val; % Table 1 left column, S2 Table left column,
                        % Fig 1B (dominating components)
                        % Fig 4A red arrows

%% Calulate the kinetics of the individual components
c_0 = [1;zeros(n_inter-1,1)]; % initial values of the concetrations
f = diag(eig_vect\c_0);
A = eig_vect*f;
c_solv = A*exp(-time./true_tau);

%% Calculate the DADSs
DADS = sp_dif(:,sp_dif_assign)*A(1:end-1,:);

%% Calculate the observable absorption kinetic data
raw_data = DADS*exp(-time./true_tau);

%% Put the data to a container
raw.raw_data = raw_data;
raw.time = time;
raw.raw_wavelength = raw_wavelength;
raw.lamp = lamp;
raw.sp_abs = sp_abs;
raw.sp_dif = sp_dif;
raw.c_solv = c_solv;
raw.DADS = DADS;

