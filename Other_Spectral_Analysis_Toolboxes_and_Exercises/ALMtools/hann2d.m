function [H Wss] = hann2d(M)

% [H Wss] = hann2d(M) 
% 
% Windows matrix M with an elliptical Hann (raised cosine) window. Returns
% the windowed data in H. Also returns the summed square of weighting 
% coefficients, Wss, used in the normalization of the power spectrum.

% (c) 2007 Taylor Perron

[ny nx] = size(M);
a = (nx-1)/2; b = (ny-1)/2; % matrix coordinates of centroid of dem
[X Y] = meshgrid(1:nx,1:ny);

warning off MATLAB:divideByZero
theta = (X==a).*(pi/2) + (X~=a).*atan2((Y-b),(X-a)); % angular polar coordinate
warning on MATLAB:divideByZero

r = sqrt((Y-b).^2 + (X-a).^2); % radial polar coordinate
rprime = sqrt((a^2)*(b^2)*(b^2*(cos(theta)).^2 + a^2*(sin(theta)).^2).^(-1)); % 'radius' of ellipse for this theta

hanncoeff = (r < rprime).*(0.5*(1 + cos(pi*r./rprime)));
H = M.*hanncoeff;

Wss = sum(sum(hanncoeff.^2));