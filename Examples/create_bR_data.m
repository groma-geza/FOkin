function [data, time, wavelength] = create_bR_data(noise_sigma,...
    wavelength_step,  show)

if nargin < 3
    show = 0;
end

raw = create_raw_bR_data;

%% Get the data from the container
raw_data = raw.raw_data;
time = raw.time;
raw_wavelength = raw.raw_wavelength;
sp_abs = raw.sp_abs;
sp_dif = raw.sp_dif;
c_solv = raw.c_solv;
DADS = raw.DADS;

%% Transpose raw_data, raw_wavelength and time for convenient use
raw_data = raw_data';
raw_wavelength = raw_wavelength';
time = time';

%% Compress the dataset by averaging over segments of wavelength
[m, p0] = size(raw_data);
p = floor(p0/wavelength_step);
assert((p > 0), 'Wavelength step is too large.')
step = p0/p;
wavelength = NaN(1, p);
data = NaN(m, p);
jj = 1;
hi = 0;
for ii = 1:step:(step*(p - 1) + 1)
    lo = hi + 1;
    hi = round(ii + step) - 1;
    % average
    wavelength(jj) = floor(mean(raw_wavelength(lo:hi)));
    data(:,jj) = mean(raw_data(:,lo:hi), 2);
    jj = jj+1;
end

%% Add Gaussian noise to the data
% Set the random number generator as it was in the presented calculation
% to ensure generating identical noise.
rng(2)
randn(1, 1900);
noise = noise_sigma*max(max(abs(data)))*randn(size(data));
data = data + noise;
rng('shuffle') % shuffle the random number generator
               % for further calculations

%% Visualize the different data
if show
    % Plot the absolute end difference spectra
    colors = [0,    0,    1;
        0,    0.5,  0;
        1,    0,    0;
        0,    0.75, 0.75;
        0.75, 0,    0.75;
        0.75, 0.75, 0;
        0.25, 0.25, 0.25];
    
    figure % S1 Fig A
    plot(raw_wavelength, sp_abs, 'LineWidth', 1)
    xlabel('Wavelength (nm)')
    ylabel('Absorbance (a. u.)')
    grid on
    legtxt = {'K1,2','L1','L2,3','M1,2,3','N1,2','O','bR'};
    l = legend(legtxt);
    l.Position(1:2) = [.29 .66];
    hax = gca;
    hax.ColorOrder = colors;
    
    figure % S1 Fig B
    plot(raw_wavelength, sp_dif, 'LineWidth', 1)
    xlabel('Wavelength (nm)')
    ylabel('Absorbance difference (a. u.)')
    grid on
    l = legend(legtxt(1:6));
    l.Position(1:2) = [.38 .66];
    hax = gca;
    hax.ColorOrder = colors;

    
    % Plot the kinetics of the individual components
    colors = [0,    0,    1;
        0,    0,    1;
        0,    0.5,  0;
        1,    0,    0;
        1,    0,    0;
        0,    0.75, 0.75;
        0,    0.75, 0.75;
        0,    0.75, 0.75;
        0.75, 0,    0.75;
        0.75, 0,    0.75;
        0.75, 0.75, 0;
        0.25, 0.25, 0.25];
    
    line_styles = {'-', '--', '-', '--', '-.', '-', '--', '-.', '-', '--',...
        '-', '-'}';
    
    figure % S1 Fig C
    hlins = semilogx(time, c_solv, 'LineWidth', 1);
    set(hlins, {'LineStyle'}, line_styles)
    xlabel('Time (s)')
    ylabel('Conentration (a. u.)')
    grid on
    leg=legend('K1','K2','L1','L2','L3','M1','M2','M3','N1','N2','O','BR');
    leg.Position(1:2) = [0.84, 0.25];
    hax = gca;
    hax.ColorOrder = colors;

    % Plot the DADSs
    figure % Fig 1A (dominating components)
    plot(raw_wavelength, DADS)
    xlabel('Wavelength')
    ylabel('Amplitude (rel)')
    grid on

    % Plot the observable absortion kinetic data in temporal and spectral
    % representation
    figure % Fig 1C
    semilogx(time, data)
    xlabel('Time (s)')
    y_label = '\DeltaA (rel)';
    ylabel(y_label)
    
    figure % Fig 1D
    plot(wavelength, data)
    xlabel('Wavelength (nm)')
    ylabel(y_label)
end
