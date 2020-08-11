clf
%Make synthetic signal from three sinusoids
x=0:0.01:15.9155*pi;
y1=sin(x);
y2=3*sin(x);
y3=sin(8*x);
y4=sin(x+1);
ytot=y1+y2+y3+y4;

%Remove 500 random indices (replace with NaN)
r=randperm(length(x),round(length(x)*0.2));
ysig=ytot;
ysig(r)=NaN;

%Calculate Lomb-Scargle periodograms
dx=.01;
fs=1/dx;
N=length(x);
[plytot,fltot]=plomb(ytot,fs);
[plysig,flsig]=plomb(ysig,fs);

py1=fft(y1);
py2=fft(y2);
py3=fft(y3);
py4=fft(y4);
pytot=fft(ytot);

fy=fs*(0:(N-2)/2)/N;


figure(1)
subplot(2,1,1)
plot(x,y1,'r',x,y2,'g',x,y3,'b',x,y4,'m',x,ytot,'k')
subplot(2,1,2)
plot(fy,py1(1:length(x)/2),'r',fy,py2(1:length(x)/2),'g',fy,py3(1:length(x)/2),'b',fy,py4(1:length(x)/2),'m',fy,pytot(1:length(x)/2),'-k',fltot,plytot,'-y',flsig,plysig,'-c')


