clear
dx=0.01;
x=0:dx:10;
sine=@(f,A,p) A.*sin(2.*pi.*f.*x+p);


y1=sine(10,1,0);
y2=sine(1/30,5,0);
y3=sine(1/5,3,0);
y4=sine(1,3,0);
y5=sine(5,5,0);

y=y1;%+y2+y3+y4+y5;


%Spectral Analysis
% Y-dimension roughness: FFT elevation y
figure(3)
clf
legtext={};
for i=8:14
    n=2^i; %number of samples to use for fft (using same number generates equal frequencies for all sites)
    Y=detrend(y); %detrend
    Y=taper(Y,0.05,0.05); %taper to avoid high frequency edge effects
    Ym=fft(Y,n).*dx; %fourier transform elevation and multiply by dy to keep real units [m]
    Pyy=2.*abs(Ym.^2); %calculate periodogram spectrum; multiply by 2 to maintain units (m), aliasing power 'folded' from higher frequencies we are ignoring
%     Pyy=Pyy./(dx*n); %normalize to power per unit frequency
    Pyy=Pyy(1:end/2); %ignore redundant frequencies (Nyquist)

    f=(0:n/2-1)./(n*dx);
    dw=min(eldiff(fliplr(1./f(2:end))));
    w=0:dw:1/f(2);
    mmwindow=0.02; %size of moving mean window in [m]
    k=mmwindow/dw;

    [Pmax,iPmax]=max(Pyy);
    fPmax=f(iPmax);
    LPmax=1/fPmax;

    Pyyint=interp1(f,Pyy,1./w); %interpolate power to smallest wavelength spacing for even moving mean windowing
    Pyymm=movmean(Pyyint,k); %moving mean to smooth spectrum

    %Plot spectra
%     subplot(1,3,1)
%     plot(x,y1,x,y,'k');%x,y2,x,y3,x,y4,x,y5,

    subplot(1,2,1)
    loglog(f,Pyy);
    hold on
    ylabel('Power (m^2)')
    xlabel('Frequency [1/m]')

    subplot(1,2,2)
    loglog(w,sqrt(Pyymm));
    hold on
    ylabel('Amplitude (m)')
    xlabel('Lengthscale of Topographic Roughness [m]')
    
    legtext=[legtext ['n = ',num2str(n),' = ',num2str(n/length(y)),' N_y']];
end
legend(legtext)
