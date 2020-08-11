function [xout,yout] = envelope(x,y)

%[xout,yout] = envelope(x,y);
%Finds the maximum y value for each value of x to determine the upper
%envelope of y plotted against x.
%
%A.M. Booth (updated 11/2008)

%Initiate output vectors:
xout = zeros(length(x),1);
yout = zeros(length(y),1);

for i = 1:length(x)
    if x(i) == 0 
        yout(i) = 0;                %gets rid of zeros
    else
        ind = find(x(i) == x);
        yout(i) = max(y(ind));      %replaces all y with max y for given x
    end
end

for i = 2:length(x)
    if x(i) == x(i-1)
        yout(i) = 0;                %keeps only one value for each x
    end
end

%Clean up the zeros to create non-redundant vectors:
ind = find(yout > 0);
xout = x(ind);
yout = yout(ind);