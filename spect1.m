function [F,T,P]=spect1(y,t,dt,varargin)


% Produces a spectrogram with the specified resolution (in units of t), where  
% the value of each spectrogram segment is its pwelch average using windows
% of length nfft with 50% overlap.
%
% INPUTS:
%   y: data array
%   t: sample  array
%   dt: sampling interval
% 
% OPTIONAL INPUTS:
%   nfft: number of discrete Fourier transform points to use in the PSD estimate. 
%         set equal to the next power of 2 greater than the sampling frequency.
%   resolution: number of seconds in each segment of final spectrogram
% 


% PARSE OPTIONAL INPUT ARGUMENTS
okargs={'nfft' 'resolution'};
defaults={[]  []};
[nfft,resolution]=internal.stats.parseArgs(okargs, defaults, varargin{:});


fs=1/dt; % sampling frequency = 1000 Hz
segment=round(fs*resolution); % resolution: 60 second windows (number of samples per 60 seconds)

% Only evaluate if entire y vector is at least one minute long
% if segment < length(y)   
%        
%     % If t does not start at the beginning of a minute, start at next one
%     tstart=datestr0(t(1),2013);
%     if isequal(tstart(18:23),'00.000')==0
%         i=1;
%         while isequal(tstart(18:23),'00.000')==0
%             i=i+1;
%             tstart=datestr0(t(i),2013);
%         end
%         t=t(i:end);
%         y=y(i:end);
%     end
%     
%     if segment < length(y)

        % Fill gaps with NaN
        [t,y]=nanfill(t,dt,y);

        % Take pwelch for each segment
        window=nfft;%(max(f)-1)*2; % pwelch window length to average over for each segment
        noverlap=fix(window/2); % overlap by 50%        
        nseg=floor(length(y)/segment); % number of segments in y, cutting off end excess
        P=zeros(nfft/2+1,nseg); % power spectrum array
        T=zeros(1,nseg); % time array
        for i=1:nseg;
            istart=((i-1)*segment)+1;
            istop=i*segment;
            [P(:,i),F]=pwelch(y(istart:istop),window,noverlap,nfft,fs);% P comes out in units of (nm/s)^2/Hz
            T(i)=t(istart-1+segment/2); %
        end
        
%     else
%         P=[];
%         T=[];
%         F=[];
%     end  
%                 
% else
%     P=[];
%     T=[];
%     F=[];
% end


%         
%             surf(T,F,log10(abs(P)),'edgecolor','none'); 
%             view(0,90);
%             axis tight;
% %             caxis([-8 12]);%[floor(min(min(log10(abs(P))))),ceil(max(max(log10(abs(P)))))]);
%             xlabel('Time (day)'); ylabel('Frequency (Hz)');
%             c=colorbar;
%             ylabel(c,'Power Spectral Density (log_{10}|P|)');
%             
            
            
            
            
            
            
            
            
            
           