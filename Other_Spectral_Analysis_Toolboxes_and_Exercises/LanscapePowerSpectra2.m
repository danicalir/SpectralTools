%% Landscape Power Spectra, version 2 done with Taylor's code, modified by me

%% If using readopentopo for smaller regions
%Z = readopentopo('filename','WV.tif','north',54.16,'south',27.6,...
%   'west',-127.7,'east',-96.05,'demtype','SRTMGL3');

%[Z,~] = reproject2utm(Z,90);
%% For larger regions, import from Arc
path_name = 'D:\ExternalHardDrive\GraduateSchool\ForearcCascadia\SRTM\SRTM_Proj.tif';
Z = GRIDobj(path_name);
clear path_name

Z = resample(Z,1000); % Resample to 1 km grid spacing

figure
Z = crop(Z,'interactive'); % IMPORTANT if using Taylor's code: Crop so that the only NaNs are ocean
Z.Z(isnan(Z.Z)) = 0; % Turn ocean NaNs to sea level. 

dx = Z.cellsize;
dy = dx;
Z = Z.Z; % Take out of TopoToolbox GRIDobj
[Ny Nx] = size(Z); % grid dimensions

% Before calculating the power spectrum, we'll make some decisions about
% the preprocessing steps we'll take. If the DEM has a non-zero mean or a
% background slope across the entire grid, this will contaminate the
% spectrum with long-wavelength signals. So our first step will be to 
% detrend the DEM by fitting a least-squares plane to the elevations and 
% then subtracting this fit.
Zo = Z; % Save the original elevations for later
%Z = Detrend(Z); % I AM STILL NOT SURE IF THIS (detrending) IS APPROPRIATE 
% FOR WHAT WE ARE DOING....

plane = Zo - Z; % Save the least-squares plane for re-trending later

% Second, the fast Fourier transform proceeds fastest if the dimensions of 
% the input matrix are integer powers of two, which we can achieve by 
% padding the DEM with zeros. 
pad = 1; % 1 means pad the data, 0 no padding.

% Third, because the edges of our DEM are not perfectly periodic, the
% spectrum can become contaminated by frequencies used to "fit" the edge
% discontinuity. We can mitigate this effect by multiplying the DEM by a
% function that tapers to zero at the edges.
window = 1; % 1 means window the data, 0 no window

%% 2D FFT

% Calculate the power spectrum using the 2D Fast Fourier Transform (FFT).
% Note that the zero padding and windowing happen inside fft2D.
[Pm fm Pv fv] = fft2D(Z,dx,dy,pad,window); 

% A few words about the output from the previous step:
% 
% We have calculated the Discrete Fourier Transform (DFT) periodogram. This
% is different from the power spectral density, which is also commonly used
% as an estimate of the power spectrum. 
%
% Pm is the 2D power spectrum matrix, and fm is the corresponding matrix of 
% radial frequencies. Pv is a vector version of the power spectrum, and fv 
% is the corresponding vector of frequencies. Units of Pm and Pv are 
% [units of Z]^2, which in our case is length^2 since Z is a matrix of 
% elevations. Units of fm and fv are [units of x]^-1, which for us is 
% 1/length since x is distance. Figs. 3a-e show examples of a simple 1D
% signal and its power spectrum, and Fig. 3f shows an example of a simple
% surface and its 2D power spectrum.
%
% What do the units mean? The periodogram is a measure of how much of the
% original elevation field's variance falls within a given frequency range.
% You can check that the sum of the periodogram is roughly equal to the
% variance in Z. (It will be somewhat less due to the zero padding.)
%
% What about the frequency limits? Wavelength is equal to the inverse of 
% frequency. The shortest wavelength that can be resolved by 2D data is 
% oriented diagonally, and is equal to sqrt(dx^2 + dy^2), or sqrt(2) for 
% our data (Fig. 4). The longest wavelength that can be resolved is 
% infinite in principle, but in practice we wouldn't feel comfortable 
% trying to detect any signal with a wavelength longer than our dataset. To 
% be even more conservative, the longest wavelength we trust should be a 
% few times shorter than the width of the DEM -- a few hundred meters in 
% this case.

% Clean up workspace
clear pad window

%% 1D Power Spectrum

% Plot the 1D version of the spectrum (Fig. 5). We'll plot the 2D spectrum 
% a few steps later, but for now it's easier to visualize in 1D.
% s1d = figure('Name','Fig. 5: 1D Spectrum','NumberTitle','off');
% subplot(2,1,1)
[axf axw] = SpecPlot1D(downsample(fv,4),downsample(Pv,4));

%% Background Spectrum

% Now let's see what the spectrum looks like in 2 dimensions. If we just
% looked at it as is, we wouldn't see much: as the 1D spectrum shows,
% longer-wavelength signals have much higher power, and would swamp any
% shorter-wavelength detail. To make the 2D spectrum easier to view, we'll
% remove the power-law background trend. Because there are so many more
% points at higher frequencies, we'll bin the 1D spectrum and fit the trend
% to the binned values (Fig. 5).

nbin = 20; % Number of bins
B = bin(log10(fv),log10(Pv),nbin,0); % Bin the log-transformed data

% Plot the binned values
hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w')

% Fit a trend with the form P ~ 1/f^n, and plot it
fit = robustfit(B(:,1),B(:,2));
plot(axf,10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k')

%% 2D Power Spectrum

fm = single(fm);  % Will help with memory on larger DEMs
% Use this fit to normalize the 2D spectrum
plots = 1;
[Pnorm Pvecnorm fvec] = fft_normpower(Pm,(10^fit(1)*fm.^fit(2)),fm,plots);
hold on
plot([1/(2*sqrt(2)*dx) 1/(2*sqrt(2)*dx)],[0 22],'r') % Nyquis wavenumber (frequency)
plot(fvec,(1.96*nanstd(Pvecnorm,0,'all'))*ones(length(fvec),1),'b'); % 95th Percentile

%% Test wavelet
% These values will need to be calculated/observed from figures above. 
peak1 = 1/(4.9*10^(-6));     % peaks in the frequency domain
peak2 = 1/(9.7*10^(-6));
peak3 = 1/(2.89*10^(-5));

a1 = (peak1/(2*pi*dx))*sqrt(5/2);   
a2 = (peak2/(2*pi*dx))*sqrt(5/2); 
a3 = (peak3/(2*pi*dx))*sqrt(5/2); % Include more values as necessary

[C1] = conv2_mexh(Zo,a1,dx); % Apply wavelet to original, non-detrended (if that was done) DEM.
[C2] = conv2_mexh(Zo,a2,dx);  
[C3] = conv2_mexh(Zo,a3,dx);

figure
subplot(2,2,1)
imagesc(Z);
title('Original Topography')
axis image

subplot(2,2,2)
imagesc(C1)
title(['Topography corresponding to ',num2str(round(peak1/100)),' km wavelength'])
axis image

subplot(2,2,3)
imagesc(C2)
title(['Topography corresponding to ',num2str(round(peak2/100)),' km wavelength'])
axis image

subplot(2,2,4)
imagesc(C3)
title(['Topography corresponding to ',num2str(round(peak3/100)),' km wavelength'])
axis image

figure
subplot(2,2,1)
imagesc(Z);
title('Original Topography')
axis image

subplot(2,2,2)
imagesc(C1.^2)
title(['Topography corresponding to ',num2str(round(peak1/100)),' km wavelength, squared'])
axis image

subplot(2,2,3)
imagesc(C2.^2)
title(['Topography corresponding to ',num2str(round(peak2/100)),' km wavelength, squared'])
axis image

subplot(2,2,4)
imagesc(C3.^2)
title(['Topography corresponding to ',num2str(round(peak3/100)),' km wavelength, squared'])
axis image
