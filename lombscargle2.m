clear
close all
clims=[];

load('Data/topo.mat');
iorder=[2 4 6 1 3 5 7]; %reorder colors from figure order from plot_rockdst_disentrainmentrate_FINAL.m
colors=[237 177 32; 179 177 32; 119 172 48; 119 172 48; 77 190 238; 0 131 205; 126 47 142]./255;%colors for plotting hi-pass
colors=colors(iorder,:);
linestyle={'--' '--' '--' '-' '-' '-' '-'};

for i = 1:7
  
    raster={topo.raster{i} topo.raster_lp{i} topo.raster_hp{i}};
   
    
    %SPECTRA
        
        %Calculate a 1D periodogram for each of these random samples
        fs = 1/.01; %units of samples/meter
        [PLz{i},fLz{i}] = plomb(raster{1}',fs);
%         [PLlp,fLlp] = plomb(raster{2}',fs);
%         [PLhp,fLhp] = plomb(raster{3}',fs);
        
        PLzm{i}=mean(PLz{i}');
%         PLlp=mean(PLlp');
%         PLhp=mean(PLhp');
        
        %turn frequency vectors into wavelength (= lengthscale of topographic roughness = 1/freq)
        wLz{i}=1./fLz{i};
%         wLlp=1./fLlp;
%         wLhp=1./fLhp;
        
        figure(1)
        subplot(2,1,1)
        loglog(wLz{i},real(PLzm{i}),'Color',colors(i,:),'LineStyle',linestyle{i},'linewidth',1);
        hold on        
        xlabel('Lengthscale of Topographic Roughness [m]')
        ylabel('PSD [m^2/m]') 
        
        subplot(2,1,2)
        loglog(fLz{i},real(PLzm{i}),'Color',colors(i,:),'LineStyle',linestyle{i},'linewidth',1);
        hold on
        xlabel('Frequency of Topographic Roughness [1/m]')
        ylabel('PSD [m^2/m]')        
    
end       
    
save('Data/plombperiodograms.mat','PLz','PLzm','fLz','wLz','Site');