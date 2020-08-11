function D = nanDetrend(M)

% D = nanDetrend(M)
%
% Fits a plane to the surface (containing NaN elements) defined by the
% elements of matrix M and detrends the surface by subtracting the value of
% the planar fit at each element. Returns the detrended surface in matrix D.
% 
% Dependencies: lsplane.m

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.
% 
% Edited by Danica Roth, July, 2020

[ny nx] = size(M);

[X Y] = meshgrid(1:nx,1:ny);

[centroid, cosines] = lsplane([X(~isnan(M)) Y(~isnan(M)) M(~isnan(M))]);

% for a plane with equation z = ax + by + c
a = -cosines(1)/cosines(3);
b = -cosines(2)/cosines(3);
c = centroid(3) + ((cosines(1)*centroid(1) + cosines(2)*centroid(2))/cosines(3));

% at each (x,y) point, subtract the value of the fitted plane from M
D = M - (a*X + b*Y + c);
