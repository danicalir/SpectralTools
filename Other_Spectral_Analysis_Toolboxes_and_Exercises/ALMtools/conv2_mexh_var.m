function [Vcwt,frq,wave] = conv2_mexh_var(patch,scales,dx)

%[Vcwt,frq,wave] = conv2_mexh_psd(patch,scales,dx);
%Computes the wavelet spectrum (wavelet variance v. scale) of
%patch over specified scales using the Mexican hat wavelet.  Note that
%edge effects increase with scale, so Vcwt at larger scales becomes
%increasingly biased towards nodes in the center of patch.  
%
%Inputs:
%patch = patch of digital elevation model
%scales = wavelet scales over which Vcwt is computed
%dx = grid spacing
%
%Outputs:
%Vcwt = vector of wavelet variance computed at each scale
%frq = vector of frequencies
%wave = vector of wavelengths
%
%A.M. Booth (updated 11/2008)

tic     %start timer

%Count number of nodes:
[nrows ncols] = size(patch);
nodes = nrows*ncols;

%Normalize patch to have unit variance:
patch = patch/std(patch(:));

%Initialize output vector:
Vcwt = zeros(1,length(scales));

%Determine extent of edge effects at largest wavelet scale sampled.  NaN
%values will be assigned to the fringe of each C grid in the loop so that
%the same number of nodes are used at each scale for determining Vcwt:
fringeval = ceil(4*max(scales));

%Start counter:
k = 0;
%Loop through scales:
for a = scales
    
    %Update counter:
    k = k + 1;                  
    
    %Compute 2D CWT by calling conv2_mexh function (below):
    C = feval(@conv2_mexh1,patch,a,dx);
    
    %Mask edge effects with NaN (no data) values.
    C(1:fringeval,:) = NaN;
    C(:,1:fringeval) = NaN;
    C(nrows-fringeval+1:nrows,:) = NaN;
    C(:,ncols-fringeval+1:ncols) = NaN;
    
    %Find NaN and replace with 0:
    ind = find(isnan(C));
    C(ind) = 0;
    
    %compute wavelet variance at current scale, using number of
    %real-valued nodes:
    Vcwt(k) = 1/(2*(nodes-length(ind)))*sum(sum((C.^2),2),1);
    
end

%frequency and wavelength vectors:
wave = 2*pi*dx*scales/(5/2)^(1/2);
frq = 1./wave;

toc     %stop timer

%-------------------------------------------------------------------------%

function [C] = conv2_mexh1(patch,a,dx)
    
%Generate the mexican hat wavelet kernel at wavelet scale a.  The kernal 
%must be large enough for the wavelet to decay to ~0 at the edges.  The 
%Mexican hat is proportional to the second derivitive of a gaussian.
[X,Y] = meshgrid(-8*a:8*a,-8*a:8*a);
psi = (1/a).*(2 - (X/a).^2 - (Y/a).^2).*exp(-((X/a).^2 + (Y/a).^2)/2);

%Convolve patch with psi using matlab's conv2 function, multiplying by dx^2
%to approximate the double integral.  'same' crops C to same size as patch.
C = (dx^2)*conv2(patch,psi,'same');