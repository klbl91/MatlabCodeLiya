%SSE fit
function SineParams=sineFitSSE(x,y,fz,varargin)%SineParams /paramsFFTp
%Purpose: Estimation of noisy sine curve parameters by FFT and non linear fitting.
% Syntax:
%       [SineParams]=sineFit(x,y,optional)
%       Input: x and y values, y=offs+amp*sin(2*pi*f*x+phi)+noise
%              optional: plot graphics if ommited. Do not plot if 0.
%       Output: SineParams(1): offset (offs)
%             
%               SineParams(2): amplitude(A)
%               SineParams(3): phi,angle shift
%               SineParams(4): frequency (f)[fixed at 10 Hz]
%               SineParams(5): MSE , if negative then SineParams are from FFT 
%       yOut=offs+amp*sin(2*pi*f*x+phi)

if nargin>3 && varargin{1}==0
  boolGraphic=false;
else
  boolGraphic=true;
end
%% FFT
pi2=2*pi;
NumSamples=length(x);
T=x(2)-x(1);
fNy=1/(2*T);%Nyquist frequency
offs=mean(y);%DC value, do not take (max(y)+min(y))/2!
y_m=y-offs;%FFT much better without offset
n = 128*2^nextpow2(NumSamples);%heavy zero padding
Y = fft(y_m,n);%Y(f)
n2=floor(n/2);
P2 = abs(Y/NumSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = (0:n2)/n/T;% frequency scale
% %FFT parameters at peak
[maxFFT,maxFFTindx]=max(P1);%Peak magnitude and location
fPeak=fs(maxFFTindx);% f at peak
Phip=angle(Y(maxFFTindx))+pi/2;%Phi-Peak is for cos, sin(90°+alpha)=cos(betta), alpha=-betta
Phip=Phip-x(1)*fPeak*pi2;%shift for phi at x=0
if Phip<0;Phip=2*pi+Phip;end

%fit function ======================================================
fitfun=fittype(@(offs,A,phi,x) offs+A*sin(2*pi*fz*x+phi) );%offs+B*x+A*sin(2*pi*10*x+phi)
[fitted_curve,gof]=fit(x,y,fitfun,'StartPoint',[offs,(max(y)-min(y))/2,Phip]);
SineParams(1:3)=coeffvalues(fitted_curve);
SineParams(4)=fz;
SineParams(5)=gof.sse;
SineParams(6)=gof.rsquare;
% if gof.rsquare<0.99
%     fprintf(2,'Warning: bad fitting result for this cycle!\nFrequency of %4.2f Hz from FFT is used\n',...
%         fPeak)
%     fitfun=fittype(@(offs,A,phi,x) offs+A*sin(2*pi*fPeak*x+phi) );%offs+B*x+A*sin(2*pi*fPeak*x+phi)
%     [fitted_curve,gof]=fit(x,y,fitfun,'StartPoint',[0,0,0]);
%     SineParams(1:3)=coeffvalues(fitted_curve);
%     SineParams(4)=fPeak;
%     SineParams(5)=gof.sse;
%     SineParams(6)=gof.rsquare;
% end

if boolGraphic
  PlotResults(x,y,SineParams);
end

%% Plot results (optional)
function PlotResults(x,y,SineParams)
xstart=x(1);
xend=x(end);
xLength=xend-xstart;
xSstep=min(xLength/100,1/10*0.1);
xS=xstart:xSstep:xend;

y5=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(4)*xS+SineParams(3));%result
hFigPlotSin = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotSin' );
if isempty(hFigPlotSin)
  screensize=get(0, 'MonitorPositions');
  hFigPlotSin=figure('Tag','Fig$PlotSin','Name','Sinus',...
    'OuterPosition',[960,screensize(1,4)/2,screensize(1,3)-960,screensize(1,4)/2]);
  drawnow
end
figure(hFigPlotSin(1));
cla reset;
plot(x,y,'k.');%time series as dots
xlabel('Time [s]');
hold on;
pIn=plot(x,y,'r-');%time series as line
pResult=plot(xS,y5,'b-');%result
legend([pIn,pResult],'Input','Result');
hold off;
grid on;

 disp(['Result:        y= ' num2str(SineParams(1)) ' + '... %num2str(SineParams(2)) ' *+' 
      num2str(SineParams(2)) ...
   ' * sin(2*pi*' num2str(SineParams(4)) '+' num2str(SineParams(3)) ')   SSE: ' num2str(SineParams(5))...
   'R2: ' num2str(SineParams(6))]);
