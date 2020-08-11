function [Pnorm Pvecnorm fvec] = fft_normpower(Phigh,Plow,fmat,plots)

%[Pnorm Pvecnorm fvec] = fft_normpower(Phigh,Plow,fmat);
%Normalizes Phigh by Plow and sorts into vectors.
%
%Phigh = 2D DFT periodogram to be normalized
%Plow = normalizing periodogram
%fmat = grid of radial frequencies sampled
%plots = 1 yes, 0 no
%
%fvec = vector from freqmat for plotting
%Pnorm = normalized 2D DFT periodogram
%Pvecnorm = normalized 1D power vector
%
%A.M. Booth (updated 04/2011)

%Normalize Phigh by Plow.
Pnorm = Phigh./Plow;

%Reshape to sorted vectors:
tpow = length(fmat);
fvec = reshape(fmat,tpow*tpow,1);
Pvecnorm = reshape(Pnorm,tpow*tpow,1);
FP = [fvec Pvecnorm];
FP = sortrows(FP,1);
fvec = FP(:,1);
Pvecnorm = FP(:,2);

%(Optional) Plot the normalized spectra:
if plots == 1
    figure; imagesc(Pnorm); axis equal tight; colorbar;
    figure; semilogx(fvec,Pvecnorm,'.k');
end