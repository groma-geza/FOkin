[data, time, wavelength] = create_bR_data(1.E-3, 10); % create test data
fokin = FOkin(data, time, wavelength); % create a FOkin object
fokin.options.signal_label % display the value of a property chain
opt =  fokin.options; % make a copy of the first property (by reference!)
opt.signal_label = '\DeltaA (rel)'; % assign a new value
fokin.options.signal_label % test it on the original object