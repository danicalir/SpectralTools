function [Vdftsum,fband] = fft_powsum(dem,w,dx,fmin,fmax,normalize,plots)

%[Vdftsum,fband] = fft_powsum(dem,w,dx,fmin,fmax,normalize);
% 
%Computes sum of spectral power (variance) at frequencis from fmin to fmax 
%in a moving window of size w x w at each point in dem.
%
%Inputs:   dem = data matrix (elevations)
%          w = width of window (must be odd)
%          dx = grid spacing
%          fmin, fmax = min and max frequency over which power is summed
%          normalize = 1 yes, 0 no
%          plots = 1yes, 0 no
%
%Outputs:  Vdftsum = Sum of DFT periodogram (variance) between fmin and 
%                    fmax at each point in dem
%          fband = characteristic frequency and wavelength band
%
%Dependencies: detrend_tp.m --> lsplane.m, hann2d.m (all below)
%
%A.M. Booth (updated 04/2011)

tic                                     %start timer

%Before starting the double loop, generate all grids and vectors needed for
%use in the loop:

[nrows ncols] = size(dem);              %dimensions of dem
center = floor(w/2)+1;                  %define center of window
Vdftsum = zeros(nrows,ncols);           %initiate output grid

%Initiate grids for detrend1.m and hann2d1.m:
[X Y] = meshgrid(1:w,1:w);

%Calculate power of 2 for zero padding:
tpow=2.^(ceil(log(w)/log(2)));

%Calculate the frequency bin size.  Frequency goes from zero (DC) to
%1/(2*dx) (Nyquist) in twopower/2 increments:
df = 1/(dx*tpow);

%Generate a frequency array for filtering the DFT periodgram in the loop:
xc = tpow/2+1; yc = xc;                 %array indices of zero frequency
[cols rows] = meshgrid(1:tpow);         %column and row indices
fmat = sqrt((df*(rows-yc)).^2 + (df*(cols-xc)).^2);

%Set entries in fmat between fmin and fmax to 1, other entries to 0, to 
%create a mask for filtering in the loop:
fmat(fmat <= fmin | fmat >= fmax) = 0;
fmat(fmat >= fmin & fmat <= fmax) = 1;

%(Optional) Filter by direction:
if 1 == 1
    %Set min and max theta from -180 to +180 degrees: 
    thetamin = 40*pi/180;               %radians
    thetamax = 48*pi/180;
    %Indices of elements in frequency matrix relative to center:
    [Xfilt Yfilt] = meshgrid(-tpow/2:1:tpow/2-1,-tpow/2:1:tpow/2-1);
    %Create grid of angles (in radians) corresponding to fmat:
    thetamat = -atan2(Yfilt,Xfilt);             %ccw is positive
    %Set entries in fmat outside theta range to zero.  Note, this does not 
    %include the redundant half of the 2D periodogram, for example 
    %specifying angles of 0 to 40 degrees does not include the redundant 
    %spectral power at angles of -140 to -180.  
    fmat(thetamat < thetamin | thetamat > thetamax) = 0;
end

%Do a 2D FFT in a moving window of size w x w:
for m = center:(nrows-center+1)
    for n = center:(ncols-center+1)
        
        %Define values in local window.        
        win = dem(m-center+1:m+center-1,n-center+1:n+center-1);
        
        %Only use patches that do not contain NaN values:
        if sum(sum(isnan(win)))==0
            
            %Remove planar trend.    
            win = feval(@detrend1,win,X,Y);          
        
            %(Optional) Normalize so data has unit variance:
            if normalize == 1
                win = win/std(win(:));
            end
            win_var = var(win(:));      %variance of detrended patch
            
            %Window locdem with Hann raised cosine window.
            [win] = feval(@hann2d1,win,X,Y,w);
            
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
            %variance of the detrended locpatch.  This corrects for the 
            %reduction in variance caused by the windowing function:        
            win = win_var*win/sum(win(:)); 
            %Filter w/ fmat:
            win = win.*fmat;
            %Sum get filtered variance and assign spectral power sum to 
            %output grid:
            Vdftsum(m,n) = sum(win(:));     
        else                            
            Vdftsum(m,n) = NaN;         %NaN if window contains NaNs
        end
            
    end    
end

%Output characteristic bands of frequencies and wavelengths:
fband = [fmin,fmax;1./fmin,1./fmax];
%Assign NaN values to fringe of output grid:
Vdftsum(1:center-1,1:ncols) = NaN; 
Vdftsum(1:nrows,1:center-1) = NaN;
Vdftsum(nrows-center+2:nrows,1:ncols) = NaN; 
Vdftsum(1:nrows,ncols-center+2:ncols) = NaN;

%(Optional): Image the spectral power sums:
if plots == 1
    figure; imagesc(log(Vdftsum)); axis equal tight; colorbar;
end

toc     %stop timer

%-------------------------------------------------------------------------%

function D = detrend1(M,X,Y)                %from J.T. Perron

[centroid,cosines] = lsplane1([X(:) Y(:) M(:)]);

% for a plane with equation z = ax + by + c
a = -cosines(1)/cosines(3);
b = -cosines(2)/cosines(3);
c = centroid(3) + ...
    ((cosines(1)*centroid(1) + cosines(2)*centroid(2))/cosines(3));

% at each (x,y) point, subtract the value of the fitted plane from dem
D = M - (a*X + b*Y + c);

%-------------------------------------------------------------------------%

function [x0,a] = lsplane1(N)               %from J.T. Perron

%calculate centroid
x0 = mean(N)';

%form matrix A of translated points
A = [(N(:, 1) - x0(1)) (N(:, 2) - x0(2)) (N(:, 3) - x0(3))];

%calculate the SVD of A
[~, S, V] = svd(A,0);

%find the smallest singular value in S and extract from V the
%corresponding right singular vector
[~,i] = min(diag(S));
a = V(:,i);

%-------------------------------------------------------------------------%

function [H] = hann2d1(M,X,Y,w)         %from J.T. Perron

%matrix coordinates of centroid of dem:
a = (w-1)/2; 
b = (w-1)/2;           

%angular polar coordinates:
theta = (X==a).*(pi/2) + (X~=a).*atan2((Y-b),(X-a)); 

%radial polar coordinates:
r = sqrt((Y-b).^2 + (X-a).^2); 
%'radius' of ellipse for this theta:
rprime = sqrt((a^2)*(b^2)*(b^2*(cos(theta)).^2 + ...
    a^2*(sin(theta)).^2).^(-1)); 

hanncoeff = (r < rprime).*(0.5*(1 + cos(pi*r./rprime)));
H = M.*hanncoeff;