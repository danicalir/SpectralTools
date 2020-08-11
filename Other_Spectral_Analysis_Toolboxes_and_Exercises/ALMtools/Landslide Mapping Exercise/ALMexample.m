%Example of how to create a deep-seated landslide map using the functions 
%in ALMtools, as described in Booth et al. (2009), Geomorphology 109,
%132-147, doi:10.1016/j.geomorph.2009.02.027

%Load the dem, set grid spacing, and create hillshade map for viewing:
load dem.mat;
dx = 1.8288;    %meters
[nrows,ncols] = size(dem);
demh = hillshade(dem,1:nrows,1:ncols);
figure; imagesc(demh); axis equal tight; colormap gray;

%%Generate power spectra representative of failed and unfailed terrain to
%%determine characteristic wavelengths of landslides in the study area.

%Select representative patches of failed and unfailed terrain:
fld = dem(223:325,227:385);  unfld = dem(370:560,180:420);

%Generate an average Fourier spectrum for each patch:
w = 65;             %window width must be odd
normalize = 0;      %1 = yes, normalizes each window to have variance of 1
plots = 1;          %1 = yes, will plot power spectrum (both 1D and 2D)
[Vdftave_fld,Vdftvec_fld,fvec,freqmat] = fft_meanspec(fld,w,dx,normalize,plots);
[Vdftave_unfld,Vdftvec_unfld] = fft_meanspec(unfld,w,dx,normalize,plots);

%Normalize the failed spectrum with the unfailed spectrum:
[Vdft_norm,Vdftvec_norm]=fft_normpower(Vdftave_fld,Vdftave_unfld,freqmat,plots);

%Plot Fourier spectra:
figure;  subplot(2,2,1);                       %failed and unfailed spectra 
loglog(fvec,Vdftvec_fld,'.r');
hold on;  loglog(fvec,Vdftvec_unfld,'.b');  xlim([.003 .5]);
title('1D Failed and Unfaliled Power Spectra');
ylabel('DFT Spectral Power (m^2)');
legend('Failed Terrain','Unfailed Terrain');
subplot(2,2,2);                                        %normalized spectrum
semilogx(fvec,Vdftvec_norm,'.k');
[fenv,Vdft_normenv] = envelope(fvec,Vdftvec_norm);    %envelope of spectrum
hold on;  semilogx(fenv,Vdft_normenv,'-k');  xlim([.003 .5]);
title('1D Normalized Power Spectrum');
ylabel('Normalized Power (DFT)');

%Redefine larger representative patches to reduce edge effects of the
%wavelet transform.  Edge effects extend ~4*waveletscale nodes:
fld = dem(223-36:325+36,227-36:385+36); %+/-36 b/c largest scale below is 9  
unfld = dem(370-36:560+36,180-36:420+36);
%Generate wavelet spectra for each patch:
scales = exp(linspace(0,2.2,20));         %20 log-spaced scales from 1 to 9
[Vcwt_fld,frq,wave] = conv2_mexh_var(fld,scales,dx);
[Vcwt_unfld] = conv2_mexh_var(unfld,scales,dx);

%Normalize the failed spectrum with the unfailed spectrum:
Vcwt_norm = Vcwt_fld./Vcwt_unfld;

%Plot wavelet spectra:
subplot(2,2,3);                                %failed and unfailed spectra
loglog(frq,Vcwt_fld,'ok');
hold on; loglog(frq,Vcwt_unfld,'dk');  xlim([.003 .5]);
xlabel('Spatial Frequency (1/m)');
ylabel('CWT Spectral Power (m^2)');
legend('Failed Terrain','Unfailed Terrain');
subplot(2,2,4);                                        %normalized spectrum 
semilogx(frq,Vcwt_norm,'o-k');  xlim([.003 .5]);
xlabel('Spatial Frequency (1/m)');
ylabel('Normalized Power (CWT)');

%%Map spatial patterns of spectral power contained in the characteristic
%%frequency band using a windowed 2D DFT algorithm, and a 2D CWT algorithm.

%Define characteristic frequency band from normalized Fourier spectrum:
fmin = 0.02;  fmax = 0.05;
%Windowed 2D DFT (takes about 5min):
[Vdftsum,fband] = fft_powsum(dem,w,dx,fmin,fmax,normalize,plots);
figure;  subplot(2,1,1); 
imagesc(log(Vdftsum)); axis equal tight; colorbar;
title('DFT Spectral Power Sum');

%Compute wavelet coefficents using the Mexican hat wavelet at scales
%corresponding to the characteristic frequency band:
[C2] = conv2_mexh(dem,2,dx);                                    %frq = 0.07
[C3] = conv2_mexh(dem,3,dx);
[C4] = conv2_mexh(dem,4,dx);
[C5] = conv2_mexh(dem,5,dx);
[C6] = conv2_mexh(dem,6,dx);                                    %frq = 0.02

%Square and sum wavelet coefficents over sampled scales at each node:
Vcwtsum = C2.^2 + C3.^2 + C4.^2 + C5.^2 + C6.^2;
subplot(2,1,2);  imagesc(log(Vcwtsum)); axis equal tight; colorbar;
title('CWT Spectral Power Sum');

%Smooth Vcwtsum to analyze broader spatial patterns:
radius = 25;
[Vcwtsumsmooth] = smooth2(Vcwtsum,radius);
subplot(2,1,2);  imagesc(log(Vcwtsumsmooth)); axis equal tight; colorbar;
title('CWT Spectral Power Sum (Smoothed)');

%%Generate landslide maps by classifying Vdftsum and Vcwtsumsmooth based on
%%a cutoff spectral power sum.  The optimal cutoff minimizes the error
%%index (Carrara et al., 1992) between the classified array and a map of
%%known landslides:

%Load the landslide map:
load lsmap.mat;

%Count correctly and incorrectly identified nodes for a series of cutoff
%spectral power sums.  Note that cutoff values are spaced logarithmically:   
bins = 200;                                                 %number of bins
[kdft,clsdft,cnsdft,ilsdft,insdft,numNaNsdft] = ...
    slide_stats(lsmap,log(Vdftsum),bins);                    %for DFT array
[kcwt,clscwt,cnscwt,ilscwt,inscwt,numNaNscwt] = ...
    slide_stats(lsmap,log(Vcwtsumsmooth),bins);              %for CWT array

%Find minimum error index:
[koptdft,errdft,optcurvedft] = ...
    slide_optcurve(kdft,lsmap,clsdft,cnsdft,ilsdft,insdft);  %for DFT array
[koptcwt,errcwt,optcurvecwt] = ...                         
    slide_optcurve(kcwt,lsmap,clscwt,cnscwt,ilscwt,inscwt);  %for CWT array

%Generate optimized landslide maps:
lsmap_dft = Vdftsum >= exp(koptdft(1));
lsmap_cwt = Vcwtsumsmooth >= exp(koptcwt(1));

%%Save arrays of interest as ASCII files to overlay in a GIS:
save dem.txt dem -ASCII;
save Vdftsum.txt Vdftsum -ASCII;
save Vcwtsum.txt Vcwtsum -ASCII;
save lsmap.txt lsmap -ASCII;
save lsmap_dft.txt lsmap_dft -ASCII;
save lsmap_cwt.txt lsmap_cwt -ASCII;