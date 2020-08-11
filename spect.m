% Produces a spectrogram with the specified resolution (in seconds), where  
% the value of each segment is its pwelch average using windows of length 
% nfft with 50% overlap.
%
% y: data array
% t: data time array
% dt: smpling interval
% nfft: number of discrete Fourier transform points to use in the PSD estimate. 
%       set equal to the next power of 2 greater than the sampling frequency.
% Resolution: number of seconds in each segment of final spectrogram
% 


function [F,T,P]=spect(y,t,dt,nfft,resolution)

fs=1/dt; % sampling frequency = 1000 Hz
segment=round(fs*resolution); % resolution: 2 second windows (number of samples per 2 seconds)
window=nfft; % pwelch window length to average over for each segment
noverlap=fix(window/2); % overlap by 50%        
nseg=floor(length(y)/segment); % number of segments in y
        
P=zeros((nfft/2)+1,nseg); % power spectrum
T=zeros(1,nseg); % time array
        
        for i=1:nseg;
            istart=((i-1)*segment)+1;
            istop=i*segment;
            [P(:,i),F]=pwelch(y(istart:istop),window,noverlap,nfft,fs);
            T(i)=t(istart-1+segment/2);
        end
                
        
            surf(T,F,log10(abs(P)),'edgecolor','none'); 
            view(0,90);
            axis tight;
%             caxis([-8 12]);%[floor(min(min(log10(abs(P))))),ceil(max(max(log10(abs(P)))))]);
            xlabel('Time (day)'); ylabel('Frequency (Hz)');
            c=colorbar;
            ylabel(c,'Power Spectral Density (log_{10}|P|)');
            %figtitle=horzcat(station,strcat(month(imonth,2),', 2011 Storm'));
            %title(figtitle);
            %h=get(0,'CurrentFigure');