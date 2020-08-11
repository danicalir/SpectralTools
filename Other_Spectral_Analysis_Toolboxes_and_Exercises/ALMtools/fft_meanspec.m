function [Pmat Pvec fvec fmat] = fft_meanspec(dem,w,dx,normalize,plots)

%[Pmat Pvec fvec fmat] = fft_meanspec(dem,w,dx,normalize,plots)
%Computes the 2D power spectrum of a moving window of size w x w at each
%point in dem, sums as it goes, then averages.  Output is average power 
%spectrum of input grid, useful for normalizing other spectra, or computing
%a representative short wavelength spectrum over a large patch of terrain.
%
%Inputs:   dem = digital elevation model
%          w = width of window (must be odd)
%          dx = grid spacing
%          normalize = 1 yes, 0 no
%Outputs:  Pmat = averaged 2D DFT periodogram
%          Pvec = averaged 1D power spectrum
%          fvec = vector of radial frequencies
%          fmat = matrix of radial frequencies
%
%Dependencies: detrend_tp.m --> lsplane.m, hann2d.m
%
%A.M. Booth (updated 04/2011)

tic                                     %start timer

[nrows ncols] = size(dem);              %dimensions of dem
%Define center of the moving window in it's own coordinates:
center = floor(w/2)+1;                      

%Calculate appropriate power of 2 for zero padding:
tpow = 2^(ceil(log(w)/log(2)));
%Initialize output grid:
Pmat = zeros(tpow,tpow);

%Calculate the frequency bin size (frequency goes from zero (DC) to
%1/(2*dx) (Nyquist) in tpow/2 increments):
df = 1/(dx*tpow);

%Create a matrix of radial frequencies:
xc = tpow/2+1;                          %array indices of zero frequency
yc = xc;                    
[cols rows] = meshgrid(1:tpow);         %column and row indices
%Frequency matrix.  Note that since fmat contains an even number of rows
%and columns, xc and yc are rounded up in row and column indices (shifted 
%down and to the right).  The first row and first column therefore contain 
%the Nyquist frequency in the x- and y-directions, while the last row and 
%column stop at one bin (df) below the Nyquist:   
fmat = sqrt((df*(rows-yc)).^2 + (df*(cols-xc)).^2);       

%Do a 2D FFT in a moving window of size w x w, summing on the go:
counter = 0;                            %start counter
for m = center:(nrows-center+1)
    for n = center:(ncols-center+1)
        
        %Define values in local window:
        win = dem(m-center+1:m+center-1,n-center+1:n+center-1);
        
        %Only use patches that do not contain NaN values:
        if sum(sum(isnan(win))) == 0
            counter = counter + 1;
        
            %Remove first order (planar) trend:
            win = detrend_tp(win);
        
            %(Optional) Normalize so data has unit variance:
            if normalize == 1
                win = win/std(win(:));
            end
            win_var = var(win(:));      %variance of detrended patch
        
            %Window with Hann raised cosine window:
            [win] = hann2d(win);
        
            %Do a 2D FFT, padding with zeros to tpow.  fftshift.m places
            %the zero frequency (DC) at the point (xc,yc) in the grid:
            win = fftshift(fft2(win,tpow,tpow));
            %Calculate the DFT periodogram [amplitude^2]:
            win = win.*conj(win)/(tpow^4);
            %Set power to zero at the zero frequency (DC).  After windowing
            %the data, its mean may no longer be zero, so this ensures that
            %the first-order trend is removed:
            win(xc,yc)=0;
            %Normalize so that the sum of the periodogram equals the 
            %variance of the detrended local patch.  This corrects for the 
            %reduction in variance caused by the windowing function:        
            win = win_var*win/sum(win(:));        
            %Sum up Pout each time through loop for averaging later:
            Pmat = Pmat + win;
        end
        
    end
end

%Divide by total number of times through loop to get mean.
Pmat = Pmat/counter;

%Generate sorted freqency and power vectors.  Note: these vectors 
%are redundant and could be reduced in size by half, but as coded below 
%they sum to the variance of the original data:  
Pvec = reshape(Pmat,tpow*tpow,1);
fvec = reshape(fmat,tpow*tpow,1);
FP = [fvec Pvec];
FP = sortrows(FP,1);
fvec = FP(:,1);
Pvec = FP(:,2);

%(Optional) Plot the 2D periodogram and 1D power spectrum:
if plots == 1
    figure; imagesc(log(Pmat)); axis equal tight; colorbar;
    figure; loglog(fvec,Pvec,'.k');
end

toc                                     %stop timer