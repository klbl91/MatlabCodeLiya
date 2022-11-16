%TestEKF(sDefType, sSGType, nLayer, sEngine, sAlgorithm)

function [A]=TestEKF(sDefType, sSGType, nLayer, sEngine, sAlgorithm)

noise = 0.00; % mm
TOL = 0.5; % relative tolerence for convergence test, absolute value is TOL*R
guessE = [1000, 1000, 1000]'; % MPa;

if nargin<2
    sSGType = 'linear';
end

if nargin<3
    nLayer = 3;
end

if nargin<4
    %sEngine = 'ELSYM5';
    sEngine = 'LEOP';
end

if nargin<5
    sAlgorithm = 'KALMAN';
end

% pavement structure
switch nLayer
case 3,
    layer = [150;300;0]; % mm
    E = [5000;250;60]; % MPa
    %E = [5000;30000;60]; % MPa
    v = 0.35;
case 4,
    layer = [152.4 203.2 254.4 0]';
    E = [3447.38 172.37 31026.40 51.71]';
    E = [5000 300 150 100]';
    guessE = [1000 200 100 80]';
    v = [0.35 0.40 0.25 0.45]';
    v = 0.35;
case 5,
    layer = [100;100;400;200;0];
    E = [10000;200;500;250;150];
    v = 0.35;
end

% deflection data type
switch upper(sDefType)
case 'MDD+RSD',
    
    fprintf('Doing MDD+RSD\n');
    
    Zout = [0;200;300;500;3000];
    bRelative = [0;1;1;1];
    
    Force = 40/2; % kN;
    pLoad = 0.72; % MPa;
    rLoad = sqrt(Force*1000/pLoad/pi); % mm
    XYload = [0 -150; 0 150];
    
case 'MDD',
    
    fprintf('Doing MDD\n');
    
    Zout = [150;300;500;3000];
    bRelative = [1;1;1];
    
    Force = 40/2; % kN;
    pLoad = 0.72; % MPa;
    rLoad = sqrt(Force*1000/pLoad/pi); % mm
    XYload = [0 -150; 0 150];
    
case 'FWD',
    
    fprintf('Testing FWD\n');
    
    rGeophone = [0 210 315 475 630 925 1535 1985];
    defData.File = 'Testing KalmanFilter Performance in FWD';
    defData.Drops = [];
    defData.Info = struct('Units','Metric','Date',datestr(now),'rGeoPhone',rGeophone);
    %        
    Force = 40; % kN;
    rLoad = 150; % mm
    pLoad = Force*1000/pi/rLoad/rLoad*1000; % kPa
    XYload = [0 0];
    Zout = 0;
    nGeophone = length(rGeophone);
    nDepth = 1;
    XYout = rGeophone(:)*[1 0]; 
    bRelative = 0;
    
    aDrop = struct('PointNo',1,'Station',0,'DropID',1,'Stress',pLoad,'rLoad',rLoad,...
                    'SurfaceTemperature',0,'AirTemperature',0,...
                    'Time',datestr(now),...
                    'Depth',0,'Def',rGeophone*0,...
                    'sDefType','FWD',...
                    'Filename','fake','DropSeq',1 ...
                    );
    %
    defData.Drops = aDrop;    
end

guessE = [ones(nLayer,1).*guessE];
XMap.idxE = [1:nLayer];
XMap.idxH = [];

switch lower(sSGType)
case 'linear',
    
    fprintf('Subgrade is Linear, i.e. n=0\n');
    nMET2 = 0;
    XMap.bHasNMET2 = 0;   
    guessNMET2 = [];
    
case 'nonlinear',
    
    if ~strcmpi(sEngine,'MET2')
        fprintf('ERROR: Please set Displacement Engine to MET2 to do nonlinear SG\n');
        return;
    end
    
    if any(strmatch(upper(sDefType),{'MDD+RSD','MDD'}))
        fprintf('ERROR: can not do nonlinear subgrade with in-depth deflection data\n');
        return;
    end    
    
    nMET2 = -0.9;
    fprintf('Subgrade is NonLinear, nonlinear exponent n = %f5.2\n', nMET2);

    XMap.bHasNMET2 = 1;    
    guessNMET2 = -0.5;

otherwise,
    fprintf('ERROR: Wrong type for Subgrade\n');
    
end

% Now get the exact solutions

tmp = KalmanBack(defData,'h',layer,'e',E,'nstep',1);
dispExact = tmp.Zhathist(:,1); % in mm
defData.Drops.Def = dispExact*1000; % in micron
nMea = 100;

dispMeas = repmat(dispExact,[1,nMea]) + normrnd(0,noise,length(dispExact),nMea);

switch upper(sAlgorithm)
case 'KALMAN',
    %A = EKF_Main('para',para,'disp',dispMeas,'guessE',guessE,'guessN',guessNMET2,'XMap',XMap,'tol',TOL,'nstep',100);
    A = KalmanBack(defData,'h',layer);

case 'GA',
    para.Z = DispEngine(paraExact);
    para.XMap = XMap;
    nLayer = length(layer);
    bounds = log(UC('mpa','psi',repmat([10,50000],[nLayer,1])));
    if XMap.bHasNMET2
        bounds = [bounds;tan([-0.8,0.8]*pi/2)];
    end
    nGenes = 300;
    nIter = 100;
    [topGenes, fitness,X] = GeneticOptimize(bounds, @CostFunctDispLET, nGenes, nIter, para);
    A.X = X';
    A.Z = para.Z*25.4;
    [err, A.Zhat] = CostFunctDispLET(topGenes(1,:),para);
    A.Zhat = A.Zhat*25.4;
end

A.E = UC('psi','mpa',exp(A.X(1:nLayer,:)));
if XMap.bHasNMET2
    A.n = atan(A.X(end,:))/(pi/2);
end

A.Zexact = dispExact;
hold on;
r = rGeophone;
zMea = A.Z(2:end);
zCal = A.Zhat(2:end);

%close all;

%plot(1./zMea./r,-r,'o',1./r./zCal,-r,'r-');

%plot(para.XYout(:,1)*25.4,reshape(dispExact,[nGeoPhone,nDepth]),'r');

fprintf('\n\n\t    Eexact   | Seed |   Eback | Pct Error | Thickness(mm)\n');
for i=1:length(layer)
    fprintf('\t%10.0f %9.0f %9.0f %12.5f %13.0f\n',E(i),A.E(i,1),A.E(i,end),100*abs((A.E(i,end)/E(i)-1)),layer(i));    
end
if XMap.bHasNMET2
    fprintf('\t%10.3f %9.3f %9.3f %12.5f %13.0f\n',nMET2,A.n(1),A.n(end),100*abs((A.n(end)/nMET2-1)), 0);    
end

A = struct('DefData',defData,'Result',A);