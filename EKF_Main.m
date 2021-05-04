function    [A]=EKF_Main(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FEM simulation for generating measurement value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extended Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input value
%Input_para=[E1,E2,E3]
fprintf('Setting up pavement structure\n');
defPara = ELSYM5_Init('LEOP','layer',[150;350;0], ...
                         'E',[5000;300;100],'v',0.35, ...
                         'nMET2',-0.0,'Zout',[0]);

%Inital estimation of Elastic modulus, and layer thickness if wanted
defGuessE = UC('MPa','psi',[100; 100; 100;]);

% Calculated displacement is transformed into measured displacements using matrix YMap
nDisp = defPara.nXYout*defPara.nZout;
YMap = eye(nDisp);

%Inital measurment error covariance
R = 0.01/25.4; % in inch

% Parse options

optionSet = {'para',    'inputPara',    defPara,                +1, [],...
             'disp',    'Z',            DispEngine(defPara),    +1, UC('mm','in'),...
             'engine',  'dispEngine',   'LEOP',                 +1, [],...
             'guesse',  'guessE',       defGuessE,              +1, UC('MPa','psi'),...
             'guessn',  'guessNMET2',   [-0.1],                 +1, [],...
             'guessh',  'guessH',       [],                     +1, UC('mm','in'),...
             'xmap',    'XMap',         [],                     +1, [],...
             'ymap',    'YMap',         YMap,                   +1, [],...
             'tol',     'TOL',          R*0.05,                 +1, R, ...
             'nstep',   'nStep',        20,                     +1, [], ...
             'bplot',   'bPlot',        1,                      +1, []};

[Options, sSpecified, unused] = ParseOptions(optionSet,varargin{:});
for i=1:length(unused)    
    fprintf('WARNING (EKF_Main): Unparsed input parameters - %s\n',num2str(unused{i}));
end

% Control informations used by MapXtoState function to map X to state variables
%   .idxE   index of E that are included in X
%   .idxH   index of h that are included in X
%   .bHasNMET2  flag for nonlinearity
%   .hRange mean and range for h when used, not needed when idxH is []

idx = strmatch('xmap',sSpecified,'exact');
if length(idx)==0
    XMap.idxE = 1:length(Options.inputPara.E);
    XMap.idxH = [];
    XMap.bHasNMET2 = Options.inputPara.nMET2~=0;
    Options.XMap = XMap;
end   

Options.inputPara.DispEngine = Options.dispEngine;

X = MapStateToX(Options);    

nInputs = length(Options.inputPara);
geoPhone = [];
for i=1:nInputs
    geoPhone = [geoPhone; Options.inputPara(i).XYout(:,1)];
end

P=0.8*eye(length(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalman Filter
fprintf('Kalman Filtering begins...\n');
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Weighted Global Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Err = 0;
eRate = 0;

if Options.bPlot
    close all;
    figure;
end

bConverge = 0;

noise = 0.01/25.4;

Xhist = X;
Z = Options.Z(:,1);
nMea = size(Options.Z,2);
Zhat = Z*0;

i = 0;
while i<Options.nStep
    i = i+1;
    
   Xhist(:,i) = X;
   
   iMea = mod(i,nMea)+1;
   Z = Options.Z(:,iMea);
   
   [X,P,Zhat(:,i)]=EKF_FEM(X,Z,Options,P,R);
   
   error_value=abs(Z-Zhat(:,end));
   Err(i) = UnitConvert('in','mm',max(error_value));
   ErrMin(i) = min(abs(Err));
   
   if i>10
       p = polyfit([1:5],ErrMin(i-4:i)/R,1);
%       p = polyfit([1:4],Err(i-3:i)/R,1);
       eRate(i) = abs(p(1));
   else
       eRate(i) = 1;
   end

   para = MapXtoState(Xhist(:,end),Options.XMap,Options.inputPara(1));

   fprintf('i=%5d\tmaxErr=%12.3E\n', i, Err(i));
   h = UC('in','mm',para.h);
   E = UC('psi','MPa',para.E);
   fprintf('\th=%12.1f mm, E=%12.1f MPa\n',[h,E]');       
   if Options.XMap.bHasNMET2
       fprintf('\tn=%12.2f\n',para.nMET2);
   end
   
   if Options.bPlot                    
       nGeoPhone = length(geoPhone);
       plot(geoPhone*25.4,reshape(Z,[nGeoPhone,length(Z)/nGeoPhone])*25.4,'o'); hold on;
       plot(geoPhone*25.4,reshape(Zhat(:,end),[nGeoPhone,size(Zhat,1)/nGeoPhone])*25.4,'-'); hold off;
       legend(num2str(UC('in','mm',para.Zout)));
       title(num2str(E'));
       drawnow
   end
   
   if max(error_value) < Options.TOL
       bConverge = 1;
       break
   else
       if eRate(i)<0.1
           bConverge = 0;
           break;
       elseif any(X>1E6)
           bConverge = 0;           
           break;
       end       
       
   end
   
   P=500*P;
   %check(i)=P(1,1);

%   pause
end

idx = i+1;

Xhist(:,idx) = X;
Zhat(:,idx)=h_z(X,Options.inputPara,Options.XMap,Options.YMap);
error_value=abs(Z-Zhat(:,end));
errPct = abs(error_value./Z)*100;

Err(idx) = UnitConvert('in','mm',max(error_value));
ErrMin(idx) = min(abs(Err));    


T_elapsed = toc;

nIter = length(Err)-1;
iMin = max(find(Err==min(Err)));

newSeq = [1:nIter];
if iMin~=nIter
    newSeq = [newSeq,iMin];
end

A.nIter = nIter;
A.X = Xhist(:,newSeq);
para = MapXtoState(Xhist,Options.XMap,Options.inputPara(1));
State.E = para.E;
State.nMET2 = para.nMET2;
State.h = para.h;
A.State = State;
A.bConverge = bConverge;
A.Err = Err(newSeq);
A.Z = UnitConvert('in','mm',Z);
A.Zhathist = UC('in','mm',Zhat);
A.Zhat = UnitConvert('in','mm',Zhat(:,iMin));
A.rGeophone = UC('in','mm',geoPhone);
%A.eRate = eRate(newSeq);

fprintf('Final Results:\n\tmaxErr=%12.3E\n', max(error_value));
fprintf('\tE=%12.0f\n',UC('psi','MPa',A.State.E(:,end)));       
if isfield(State,'nMET2')
    fprintf('\tn=%12.2f\n',A.State.nMET2(end));
end
fprintf('Total Number of Iterations: %5d\n',nIter);
fprintf('Total Time Used: %10.2f Seconds\n', T_elapsed);
       
