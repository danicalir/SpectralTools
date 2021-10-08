% Plots high-pass and low-pass roughness distributions and slope dependence
% for burned and vegetated sites in the OCR.

clear
close all

load('Data/topo');

iorder=[2 4 6 1 3 5 7]; %reorder colors from figure order from plot_rockdst_disentrainmentrate_FINAL.m
colors=[237 177 32; 179 177 32; 119 172 48; 119 172 48; 77 190 238; 0 131 205; 126 47 142]./255;%colors for plotting hi-pass
colors=colors(iorder,:);

linestyle={'--' '--' '--' '-' '-' '-' '-'};

iburn=1:3; %indices for burned slopes
iveg=4:7; %indices for vegetated slopes

hplegtext={};
lplegtext={};

% NON-INTERPOLATED RASTER (1 cm) 
for i=1:length(linestyle)
    [topo.CDF_nlp{i},topo.nlp{i}]=ecdf(2.*abs(topo.raster_lp{i}(:)));
    topo.n50(i)=prctile(2.*abs(topo.raster_lp{i}),50,'all'); %find median macroroughness height, n50

    figure(6)
    shp=subplot(1,2,1);
    ahp(i)=semilogx(topo.dhp{i},topo.CDF_dhp{i});
    set(ahp(i),'Color',colors(i,:),'LineStyle',linestyle{i});
    title(['High-Pass Micro-Roughness Distributions (Gaussian Filter \sigma = ',num2str(filtsigma/100),' m)']);
    ylabel('Cumulative Fraction');
    xlabel('Micro-Roughness Height d (m)');
    hold on
    
    slp=subplot(1,2,2);
    alp(i)=semilogx(topo.nlp{i},topo.CDF_nlp{i});
    set(alp(i),'Color',colors(i,:),'LineStyle',linestyle{i});
    title(['Low-Pass Macro-Roughness Distributions (Gaussian Filter \sigma = ',num2str(filtsigma/100),' m)']);
    xlabel('Macro-Roughness Height (m)');
    hold on
    
    %Text for legends to be added later
    hplegtext=[hplegtext ['S=',num2str(topo.slope(i),2),', d_{50,hp}=',num2str(topo.d50(i,1),2)]];
    lplegtext=[lplegtext ['S=',num2str(topo.slope(i),2),', n_{50,lp}=',num2str(topo.n50(i),2)]];
end
figure(6)
subplot(1,2,1)
lhp=legend(hplegtext,'Location','NorthWest');
subplot(1,2,2)
llp=legend(lplegtext,'Location','NorthWest');
set(lhp,'Box','off')
set(llp,'Box','off')
set(shp,'xlim',[0.0001 1])
set(slp,'xlim',[0.0001 10])
set(alp,'LineWidth',1)
set(ahp,'LineWidth',1)


figure(7)
clf
subplot(1,2,1)
hhpv=scatter(topo.slope(iveg),topo.d50(iveg),75,colors(iveg,:),'^','filled','SizeData',200);
hold on
hhpb=scatter(topo.slope(iburn),topo.d50(iburn),75,colors(iburn,:),'o','SizeData',200);
xlabel('slope')
ylabel('Micro-Roughness Height d [m]');

subplot(1,2,2)
lhpv=scatter(topo.slope(iveg),topo.n50(iveg),150,colors(iveg,:),'^','filled','SizeData',200);
hold on
lhpb=scatter(topo.slope(iburn),topo.n50(iburn),150,colors(iburn,:),'o','SizeData',200);
xlabel('slope')
ylabel('Macro-Roughness Height [m]');
