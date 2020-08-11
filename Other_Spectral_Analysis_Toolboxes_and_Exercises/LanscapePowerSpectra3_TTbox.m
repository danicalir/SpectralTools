%% Landscape Power Spectra and Wavelet Transform
% Compiled from code written by Taylor Perron and Adam Booth
% By Will Struble, last modified 06/28/2019
%
%% If using readopentopo for smaller regions, resolution as fine as 30m
 Z = readopentopo('filename','WV.tif','north',44.14,'south',43.67,...
   'west',-123.68,'east',-122.87,'demtype','SRTMGL1');
% 
[Z,~] = reproject2utm(Z,30);
%% For larger regions, import from Arc
path_name = 'D:\ExternalHardDrive\GraduateSchool\ForearcCascadia\SRTM\SRTM2_Proj.tif'; % DEM file
Z = GRIDobj(path_name);
clear path_name

%% Resample and crop if necessary. REMOVE ALL NaNs!
Z = resample(Z,1000); % Resample to ~1 km grid spacing for continent scale analysis

figure
Z = crop(Z,'interactive'); % IMPORTANT if using Taylor's code: Crop so that the only NaNs are ocean
Z.Z(isnan(Z.Z)) = 0; % Turn ocean NaNs to sea level. 

dx = Z.cellsize;
dy = dx;
% Z = Z.Z; % Take out of TopoToolbox GRIDobj
[Ny Nx] = size(Z.Z); % grid dimensions


Zo = Z; % Save the original elevations for later
Z = Z.Z;

% The fast Fourier transform proceeds fastest if the dimensions of 
% the input matrix are integer powers of two, which we can achieve by 
% padding the DEM with zeros. 
pad = 1; % 1 means pad the data, 0 no padding.

% Because the edges of our DEM are not perfectly periodic, the
% spectrum can become contaminated by frequencies used to "fit" the edge
% discontinuity. We can mitigate this effect by multiplying the DEM by a
% function that tapers to zero at the edges.
window = 1; % 1 means window the data, 0 no window

%% 2D FFT

% Calculate the power spectrum using the 2D Fast Fourier Transform (FFT).
% Note that the zero padding and windowing happen inside fft2D.

[Pm fm Pv fv] = fft2D(Z,dx,dy,pad,window); 

% Clean up workspace
clear pad window

%% 1D Power Spectrum
figure % SpecPlot1D does not automatically open a figure and will throw an error if the formatting is incorrect. 
[axf axw] = SpecPlot1D(downsample(fv,4),downsample(Pv,4)); %Recommend downsampling if working with large dataset
%[axf axw] = SpecPlot1D(fv,Pv);
%% Background Spectrum

nbin = 25; % Number of bins
B = bin(log10(double(fv)),log10(double(Pv)),nbin,0); % Bin the log-transformed data

% Plot the binned values
hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w')

% Fit a trend with the form P ~ 1/f^n, and plot it
fit = robustfit(B(:,1),B(:,2));
plot(axf,10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k')

Pmn = Pm./(10^fit(1)*fm.^fit(2));

%% Alternative way to normalize - Will likely lose important peaks in the spectrum though. May isolate small peaks. 
% P = polyfit(10.^B(:,1),10.^B(:,2),13);
% f = interp1(10.^B(:,1),10.^B(:,2),fv,'pchip','extrap'); % Interpolates using a 
% % shape-preserving piecewise cubic interpolation
% plot(axf,fv,f,'k.')
% 
% Pmn = Pv./f;
% figure
% semilogx(fv,Pmn,'.k');
%% 2D Power Spectrum

%fm = single(fm);  % Will help with memory on larger DEMs
% Use this fit to normalize the 2D spectrum
plots = 1;
[Pnorm Pvecnorm fvec] = fft_normpower(Pm,(10^fit(1)*fm.^fit(2)),fm,plots);
hold on
plot([1/(2*sqrt(2)*dx) 1/(2*sqrt(2)*dx)],[0 22],'r') % Nyquist wavenumber (frequency)
plot(fv,(1.96*nanstd(Pmn,0,'all'))*ones(length(fv),1),'b'); % 95th Percentile; This may not be the appropriate way

axf = gca;
frange = get(axf,'xlim');
wrange = 1./fliplr(frange);

axw=axes('Position',get(axf,'Position'),...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none',...
         'XColor','k','YColor','k',...
         'xlim',wrange,'xdir','reverse',...
         'ytick',[],'xscale','log');


set(axf,'box','off','tickdir','out')
set(axw,'box','off','tickdir','out')

set(get(axw,'xlabel'),'string','Wavelength (m)')
set(get(axf,'xlabel'),'string','Radial frequency (m^{-1})')
set(get(axf,'ylabel'),'string','Normalized Spectral Amplitude')

%% PAUSE HERE AND FIND SPECTRAL PEAKS

%% Wavelets
% These values will need to be calculated/observed from figures above. 
peak1 = 1/(4.8342*10^(-5));     % peaks in the wavenumber domain
peak2 = 1/(2.4937*10^(-5));
peak3 = 1/(5.0331*10^(-6));

%% Mexican Hay Wavelet to visualize areas of convex & convave topography

% Define scale of wavelet
a1mex = (peak1/(2*pi*dx))*sqrt(5/2);   % Torrence and Compo (1998)
a2mex = (peak2/(2*pi*dx))*sqrt(5/2); 
a3mex = (peak3/(2*pi*dx))*sqrt(5/2); % Include more values as necessary

C1mex = Zo;   % Make values compatibile w/ TopoToolbox
C2mex = Zo;
C3mex = Zo;

[C1mex.Z] = conv2_mexh(Zo.Z,a1mex,dx); % Apply wavelet to original DEM.
[C2mex.Z] = conv2_mexh(Zo.Z,a2mex,dx);  
[C3mex.Z] = conv2_mexh(Zo.Z,a3mex,dx);

%% Smooth the topography by convolving it with a Gaussian "wavelet" (not a true wavelet)
%  over the same relevant wavelength

% Define scale of Gaussian (differs from Mexican Hat for same wavelength)
a1gaus = (peak1/(2*pi*dx))*sqrt(1/2);   % Torrence and Compo (1998)
a2gaus = (peak2/(2*pi*dx))*sqrt(1/2); 
a3gaus = (peak3/(2*pi*dx))*sqrt(1/2); % Include more values as necessary

C1gaus = Zo;   % Make values compatibile w/ TopoToolbox
C2gaus = Zo;
C3gaus = Zo;

[C1gaus.Z] = conv2_gaus(Zo.Z,a1gaus,dx); % Smooth topography with Gaussian.
[C2gaus.Z] = conv2_gaus(Zo.Z,a2gaus,dx);  
[C3gaus.Z] = conv2_gaus(Zo.Z,a3gaus,dx);

%% Plot results
figure
set(gcf,'Position',get(0,'Screensize'));
subplot(1,2,1)
imagesc(C1mex);
title(['Mexican Hat Wavelet corresponding to ',num2str(round(peak1/1000)),' km wavelength'])

subplot(1,2,2)
imagesc(C1gaus);
title(['Gaussian smoothed topography, ',num2str(round(peak1/1000)),' km wavelength'])

figure
set(gcf,'Position',get(0,'Screensize'));
subplot(1,2,1)
imagesc(C2mex);
title(['Mexican Hat Wavelet corresponding to ',num2str(round(peak2/1000)),' km wavelength'])

subplot(1,2,2)
imagesc(C2gaus);
title(['Gaussian smoothed topography, ',num2str(round(peak2/1000)),' km wavelength'])

figure
set(gcf,'Position',get(0,'Screensize'));
subplot(1,2,1)
imagesc(C3mex);
title(['Mexican Hat Wavelet corresponding to ',num2str(round(peak3/1000)),' km wavelength'])

subplot(1,2,2)
imagesc(C3gaus);
title(['Gaussian smoothed topography, ',num2str(round(peak3/1000)),' km wavelength'])

%% Route Flow across topography: this will likely onlt be done for long wavelength transforms
% Modern, unfiltered topography. Be cognizant of internally drained basins. 

% Remove values of 0, to prevent flow from routing way out into the ocean.
C1plot = C1gaus;
C2plot = C2gaus;
C3plot = C3gaus;

C1plot.Z(C1plot.Z==0) = NaN;
C2plot.Z(C2plot.Z==0) = NaN;
C3plot.Z(C3plot.Z==0) = NaN;

Zoplot = Zo;
Zoplot.Z(Zoplot.Z==0) = NaN;

FD = FLOWobj(Zoplot,'preprocess','carve');
A = flowacc(FD);
S = STREAMobj(FD,A>5000);

FD1_noc = FLOWobj(C1plot,'preprocess','none');
A1_noc = flowacc(FD1_noc);
S1_noc = STREAMobj(FD1_noc,A1_noc>5000);

FD1c = FLOWobj(C1plot,'preprocess','carve');
A1c = flowacc(FD1c);
S1c = STREAMobj(FD1c,A1c>5000);

FD2_noc = FLOWobj(C2plot,'preprocess','none');
A2_noc = flowacc(FD2_noc);
S2_noc = STREAMobj(FD2_noc,A2_noc>5000);

FD2c = FLOWobj(C2plot,'preprocess','carve');
A2c = flowacc(FD2c);
S2c = STREAMobj(FD2c,A2c>5000);

FD3_noc = FLOWobj(C3plot,'preprocess','none');
A3_noc = flowacc(FD3_noc);
S3_noc = STREAMobj(FD3_noc,A3_noc>5000);

FD3c = FLOWobj(C3plot,'preprocess','carve');
A3c = flowacc(FD3c);
S3c = STREAMobj(FD3c,A3c>5000);

%% Visualize River Networks
figure
set(gcf,'Position',get(0,'Screensize'));
subplot(3,2,1)
imagesc(C1gaus)
hold on
plot(S1_noc,'k');
title(['Non-carved rivers corresponding to ',num2str(round(peak1/1000)),' km wavelength']);

subplot(3,2,2)
imagesc(C1gaus)
hold on
plot(S1c,'k');
title(['Carved rivers corresponding to ',num2str(round(peak1/1000)),' km wavelength']);

subplot(3,2,3)
imagesc(C2gaus)
hold on
plot(S2_noc,'k');
title(['Non-carved rivers corresponding to ',num2str(round(peak2/1000)),' km wavelength']);

subplot(3,2,4)
imagesc(C2gaus)
hold on
plot(S2c,'k');
title(['Carved rivers corresponding to ',num2str(round(peak2/1000)),' km wavelength']);

subplot(3,2,5)
imagesc(C3gaus)
hold on
plot(S3_noc,'k');
title(['Non-carved rivers corresponding to ',num2str(round(peak3/1000)),' km wavelength']);

subplot(3,2,6)
imagesc(C3gaus)
hold on
plot(S3c,'k');
title(['Carved rivers corresponding to ',num2str(round(peak3/1000)),' km wavelength']);

%% Extra Plotting Options

% figure
% set(gcf,'Position',get(0,'Screensize'));
% subplot(2,2,1)
% imagesc(Z);
% title('Original Topography')
% axis image
% 
% subplot(2,2,2)
% imagesc(C1)
% title(['Topography corresponding to ',num2str(round(peak1/1000)),' km wavelength'])
% axis image
% 
% subplot(2,2,3)
% imagesc(C2)
% title(['Topography corresponding to ',num2str(round(peak2/1000)),' km wavelength'])
% axis image
% colormap(cmapping2);
% colorbar
% colormap(flipud(cmapping2));
% hold on
% plot(test(:,1),test(:,2),'k','LineWidth',2)
% 
% subplot(2,2,4)
% imagesc(C3)
% title(['Topography corresponding to ',num2str(round(peak3/1000)),' km wavelength'])
% axis image
% colormap(cmapping2);
% colorbar
% colormap(flipud(cmapping2));
% 
% figure
% set(gcf,'Position',get(0,'Screensize'));
% subplot(2,2,1)
% imagesc(Z);
% title('Original Topography')
% axis image
% 
% subplot(2,2,2)
% imagesc(C1.^2)
% title(['Topography corresponding to ',num2str(round(peak1/1000)),' km wavelength, squared'])
% axis image
% 
% subplot(2,2,3)
% imagesc(C2.^2)
% title(['Topography corresponding to ',num2str(round(peak2/1000)),' km wavelength, squared'])
% axis image
% 
% subplot(2,2,4)
% imagesc(C3.^2)
% title(['Topography corresponding to ',num2str(round(peak3/1000)),' km wavelength, squared'])
% axis image
% 
% %% In case you want to see the plots on top of topography:
% figure
% subplot(2,2,1)
% imageschs(Zo);
% title('Original Topography')
% axis image
% 
% subplot(2,2,2)
% imageschs(Zo,C1)
% title(['Topography corresponding to ',num2str(round(peak1/1000)),' km wavelength'])
% axis image
% 
% subplot(2,2,3)
% imageschs(Zo,C2)
% title(['Topography corresponding to ',num2str(round(peak2/1000)),' km wavelength'])
% axis image
% 
% subplot(2,2,4)
% imageschs(Zo,C3)
% title(['Topography corresponding to ',num2str(round(peak3/1000)),' km wavelength'])
% axis image