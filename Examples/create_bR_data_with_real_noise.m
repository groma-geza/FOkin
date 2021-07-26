function [data, time, wavelength, aver_noise_level] = ...
    create_bR_data_with_real_noise(wavelength_step,  show)

if nargin < 2
    show = 0;
end

raw = create_raw_bR_data;

%% Get the data from the container
raw_data = raw.raw_data;
time = raw.time;
raw_wavelength = raw.raw_wavelength;
lamp = raw.lamp;

%% Generate the spectral distribution of the noise
sqrt_rec_lamp=1./sqrt(lamp/max(lamp)); % reciprocal of the square root of
                                       % the relative intensity spectrum of
                                       % the lamp.

% Set the random number generator as it was in the presented calculation
% to ensure generating identical noise.
rng(2)
randn(1, 1900);

[n,m]=size(raw_data);
noise = random('Poisson',1000,n,m)/1000-1; % photon noise of Poisson 
                                           % distribution

rng('shuffle') % shuffle the random number generator
               % for further calculations

scaled_noise=diag(sqrt_rec_lamp)*noise; 

%% Generate the temporal noise segments
t1=length(find(time<2e-6));
dt1 = 0.1; % dwell time in microseconds
n1=50; % number of repetitions
t2=length(find(time<2e-5)); 
dt2 = 1;
n2=25; 
t3=length(find(time<1e-4));
dt3 = 5;
n3=15; 
dt4 = 20;
n4=10; 

f1 = dt1*n1; f2 = dt2*n2; f3 = dt3*n3; f4 = dt4*n4;

scaled_noise_orig = scaled_noise;
scaled_noise(:,t1+1:t2)=scaled_noise(:,t1+1:t2)/sqrt(f2/f1);
scaled_noise(:,t2+1:t3)=scaled_noise(:,t2+1:t3)/sqrt(f3/f1);
scaled_noise(:,t3+1:end)=scaled_noise(:,t3+1:end)/sqrt(f4/f1);

factor = 0.2; % empirical factor to obtain realistic level of noise     
raw_factored_noise = factor*scaled_noise;
raw_factored_noise_orig = factor*scaled_noise_orig;

%% Transpose variables for convenient use
raw_data = raw_data';
raw_wavelength = raw_wavelength';
time = time';
raw_factored_noise = raw_factored_noise';
raw_factored_noise_orig = raw_factored_noise_orig';

%% Visualize the data
if show
    figure % S2 Fig A
    yyaxis left
    plot(raw_wavelength,lamp) 
    xlabel('Wavelength (nm)')
    ylabel('Lamp intensity (a. u.)')
    yyaxis right
    plot(raw_wavelength, sqrt_rec_lamp)
    ylabel('Noise scale factor (a. u.)')
    
    % Pick the noise temporal distribution in the absorpion peak of bR.
    wv_sel = 550;
    noise_temp_orig = raw_factored_noise_orig(:,raw_wavelength == wv_sel);
    noise_temp = raw_factored_noise(:,raw_wavelength == wv_sel);
    
    figure % S2 Fig B
    semilogx(time, noise_temp_orig, '.-b')
    hold on
    semilogx(time(t1:t2), noise_temp(t1:t2), '.-g')
    semilogx(time(t2:t3), noise_temp(t2:t3), '.-c')
    semilogx(time(t3:end), noise_temp(t3:end), '.-r')
    xlabel('Time (s)')
    ylabel('Noise')
    title(['Noise at ', num2str(wv_sel), ' nm'])
    legend('original', 'segment #2', 'segment #3', 'segment #4')    
end

%% Compress the dataset and noise by averaging over segments of wavelength
[m, p0] = size(raw_data);
p = floor(p0/wavelength_step);
assert((p > 0), 'Wavelength step is too large.')
step = p0/p;
wavelength = NaN(1, p);
data = NaN(m, p);
noise = NaN(m, p);
jj = 1;
hi = 0;
for ii = 1:step:(step*(p - 1) + 1)
    lo = hi + 1;
    hi = round(ii + step) - 1;
    % average
    wavelength(jj) = floor(mean(raw_wavelength(lo:hi)));
    data(:,jj) = mean(raw_data(:,lo:hi), 2);
    noise(:,jj) = mean(raw_factored_noise(:,lo:hi), 2);
    jj = jj+1;
end

data = data + noise;
aver_noise_level = std(noise, 0, 'all');