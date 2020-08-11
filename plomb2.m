function [Pxy,Fxy,pxy,fxy] = plomb2(Z,varargin)
%FFT2 Two-dimensional discrete Fourier Transform.
%   FFT2(X) returns the two-dimensional Fourier transform of matrix X.
%   If X is a vector, the result will have the same orientation.
%
%   FFT2(X,MROWS,NCOLS) pads matrix X with zeros to size MROWS-by-NCOLS
%   before transforming.
%

[py,fy] = plomb(Z,varargin{:});
[Pxy,fx] = plomb(py.',varargin{:});
Pxy = Pxy.';
fx = fx.';

[Fx Fy] = meshgrid(fx,fy); % matrices of column and row indices
Fxy = sqrt((Fx).^2 + (Fy).^2); % frequency matrix

fpsort = sortrows(horzcat(Fxy(:),Pxy(:)),1);
fxy = fpsort(:,1);
pxy = fpsort(:,2);