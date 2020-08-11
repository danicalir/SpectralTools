% Plots a spectrogram from input time (T), frequency (F), and power spectral density (P) arrays.
% figtitle is string (ie, 'EB1.FHE.NTE237-238')

function [h] = spectplot(T,F,P,figtitle)
surf(T,F,log10(abs(P)),'edgecolor','none');
view(0,90);
axis tight;
xlabel('Time'); ylabel('Frequency (Hz)');
c=colorbar('SouthOutside');
xlabel(c,'Power Spectral Density (log_{10}|P|)','Interpreter','Tex');
if exist('figtitle')==1
    title(figtitle);
end
h=get(0,'CurrentFigure'); 