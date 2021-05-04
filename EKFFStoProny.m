%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 05/03/2021 LJ
% Goal: use this EKF-WGI function to fit the Generalized Maxwell Model (6
% branches with one elastic single spring and 5 maxwell models) with the
% frequency sweep rest results.
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read the FS tests data
cd 'C:\1.study\PHD\Thesis\Dissertation\Chapter7-thermal cracking\Aging';
testdatapath="C:\1.study\PHD\UCPRC\Project\4.64\FAMtestresult\upscaletoHMA\fs\HMA\4PB Freq Sweep";
[basename,folder]=uigetfile(testdatapath+"\*.xlsx",'select the input data file');
fullFileName=fullfile(folder,basename);
[data,text]=xlsread(fullFileName,'Data.Clean');
Freq=data(:,1);
Stiff=data(:,4);
Phaseangle=data(:,5);
Temp=data(:,6);
StorageStiff=Stiff.*cos(Phaseangle/180*pi());
Lossstiff=Stiff.*sin(Phaseangle/180*pi());

%% EKF-WGI process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Extended Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input value
%Input_para=[E1,E2,E3,E4,E5,Einf,alfaT]
%Initial estimation of modulus
defGuessPara=[1000;1000;1000;1000;1000;1000;2];

%calculated complex modulus is transformed to measured complex modulus
%using matrix YMap
numMeasure=length(Stiff);
YMap=eye(numMeasure);

%Initial measurement error covariance
R=1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kalman Filter
fprintf('   ====================   Kalman filtering begins...   ========================\n');
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Weighted Global Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Err = 0;
eRate = 0;
bConverge=0;
noise=1e-4;
XMap.idxE=1:lengt(

%% calculate stiffness
function [Estore,Eloss]=Calstiff(Data,GuessPara)
  lamda=[0.0005,0.005,0.05,0.5,5];
  W=Data(:,1)*2*pi();
  Tr=Data(:,6)*9/5+491.67;%Rankline temperature
  Wred=10^(log10(W)+GuessPara(7)*(10^(10.5254-3.5*log10(Tr))-10^(10.5254-3.5*...
      log10(20*9/5+491.67))));
  Estore=GuessPara(6);
  Eloss=0;
  for i=1:5
    Ei1=GuessPara(i)*(Wred*lamda(i)).^2./(1+(Wred*lamda(i)).^2);
    Ei2=GuessPara(i)*(Wred*lamda(i))./(1+(Wred*lamda(i)).^2);
    Estore=Estore+Ei1;
    Eloss=Eloss+Ei2;
  end
 
end
