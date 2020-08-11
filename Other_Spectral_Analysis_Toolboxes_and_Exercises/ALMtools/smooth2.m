function [Msmooth,Mres] = smooth2(M,radius) 

%[Msmooth,Mres] = smooth2(M,radius);
%Simple smoothing function that averages array M using all points
%within specified radius.
%
%M = array to be smoothed
%radius = smoothing radius
%
%Msmooth = smoothed array M
%Mres = residuals (M - Msmooth)
%
%A.M. Booth (updated 11/2008)

tic     %start timer

%Define circular kernel to use for smoothing:  
kernel = zeros(2*radius+1,2*radius+1);
[X,Y] = meshgrid(-radius:radius,-radius:radius);
kernel((X.^2 + Y.^2).^(1/2) <= radius) = 1;

%Filter M with kernel using conv2 function. 
Msmooth = conv2(M,kernel,'same');
%Divide by the sum of all nodes in kernal to get an average at each node:
Msmooth = Msmooth/(sum(kernel(:)));

%Subtract to get residuals:
Mres = M - Msmooth;

toc     %stop timer