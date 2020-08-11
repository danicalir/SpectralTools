function [k,cls,cns,ils,ins,numNaNs] = slide_stats(mapped,identified,bins,low,high)

%[k,cls,cns,ils,ins,numNaNs] =
%slide_stats(mapped,identified,bins,low,high);
%Steps through series of cutoff values to classify identified, and counts
%correct and incorrect nodes as compared to mapped.  
%
%mapped = map of 1's and 0's indicating locations of known landslides
%identified = array of spectral power sums produced with 2D DFT or 2D CWT
%bins = number of cutoff value intervals
%low = lowest cutoff value
%high = highest cutoff value
%
%k = vector of cutoff values
%cls = correctly identified landslide nodes
%cns = correctly identified non-landslide nodes
%ils = incorrectly identified landslide nodes
%ins = incorrectly identified non-landslide nodes
%numNaNs = number of NaNs in array
%
%A.M. Booth (updated 11/2008)

tic     %start timer

%Size of array:
[nrows ncols] = size(identified);

%Find NaNs in either array, assign to mapped.  NaNs not counted:
mapped(isnan(identified)|isnan(mapped)) = NaN;
numNaNs = sum(sum(isnan(mapped),2),1);

%Check min and max, determine step size:
if nargin <=3
    disp('Min and max values in identified will be used as low and high')
    %Find minimum and maximum values in identified
    low = min(min(identified));
    high = max(max(identified));
end

%Step size:
h = (high-low)/bins;

%Initialize slidemap and output vectors:
slidemap = zeros(nrows,ncols);
cls = zeros(1,bins+1);
cns = zeros(1,bins+1);
ils = zeros(1,bins+1);
ins = zeros(1,bins+1);

%Start counting ind in output vectors:
ind = 0;

for k = low:h:high
    
    ind = ind + 1;  %count times through loop
    
    %Set all values >= k to 1, all values < k to 0:
    slidemap = identified >= k;
    
    cls(ind) = sum(sum(slidemap == 1 & mapped == 1,2),1);
    cns(ind) = sum(sum(slidemap == 0 & mapped == 0,2),1);
    ils(ind) = sum(sum(slidemap > mapped,2),1);
    ins(ind) = sum(sum(slidemap < mapped,2),1);
    
end
k = (low:h:high)';

toc     %stop timer