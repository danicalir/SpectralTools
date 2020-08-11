clear
[Z,R] = arcgridread('C:\Users\danicar\Desktop\Analysis\TopographicData\OCR\NOBLE1\Polydata_all_oct0.01m_noveg_oct0.01m_UTMz10NAD83_rotate-66_0.01krig_0.01radius_box.asc');

dx=0.01; %sampling interval 

%%Subfunction 1
for i = 1:size(Z,2)
    xfl(i,:) = [find(isnan(Z(:,i))==0,1,'first'),find(isnan(Z(:,i))==0,1,'last')]; %first and last non-NaN samples
    deltax(i) = (xfl(i,2)-xfl(i,1))*0.01; %spatial distance between first and last non-NaN values [m]
    N(i) = length(find(isnan(Z(:,i))==0)); %number of non-NaN samples [#]
end
xs = deltax./(N-1); %average sample distance [m]
fmin = 1./(4.*N.*xs); %minimum resolvable frequency [cycles/m]
fmax = 1./(2.*xs); %maximum resolvable frequency [cycles/m]
fmax_range=[min(fmax) max(fmax)];
fmin_range=[min(fmin) max(fmin)];
% Lmin = 1./fmax;
% Lmax = 1./fmin;
% Lmax_range=[min(Lmax) max(Lmax)]
% Lmin_range=[min(Lmin) max(Lmin)]

fvec(i,:)=min(fmin)


for i = 1:size(Z,2)
    [pxx{i},fx{i}] = plomb(Z(i,:),1/dx);
    [pyy{i},fy{i}] = plomb(Z(:,i),1/dx);
 
    if i==1
        L=length(fi);
    end
    
    if length(fi) <L
        
    end
end


hold on
for i=1:length(pxx)
    scatter(ones(1,length(fx{i}))*i,fx{i},8,pxx{i});
end




