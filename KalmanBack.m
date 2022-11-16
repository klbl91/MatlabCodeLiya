function [A,Eo] = KalmanBack(A,varargin)

if nargin==0    
    A = struct('File','ERROR: I don'' see any data!!!','Drops',[],'Info',[]);
    A.Info.rGeoPhone = 0;
    KalmanBack(A,'help');
    return;
end

fprintf('%s\n',A.File);

Drops = A.Drops;
Info = A.Info;
defIdxGeoPhone = 1:length(Info.rGeoPhone);
clear A;

% Default parameters

% Options allowed:

optionSet = {'nstep',       'nStep',        100,                 +1, [],...
             'n',           'guessNMET2',   0,                   +1, [],...
             'idrop',       'iDrop',        1:length(Drops),     +1, [],...
             'ignore',      'iIgnore',      [],                  +1, [],...
             'resample',    'bResample',    0,                   +0, 1,...
             'dsample',     'dSample',      100,                 +1, [],...
             'h',           'h',            [125;349;0],         +1, [],...
             'e',           'E',            [1000;200;80],       +1, [],...
             'xyload',      'XYLoad',       [0,0],               +1, [],...
             'fixede',      'iFixedE',      [],                  +1, [],...
             'fixed',       'iFixedE',      [],                  +1, [],...
             'varh',        'iVarH',        [],                  +1, [],...
             'erange',      'eRange'        [],                  +1, [],...
             'linear'       'bLinear',      1,                   +0, 1,...
             'leop',        'dispEngine',   'LEOP',              +0, 'LEOP',...
             'leopactivex', 'dispEngine',   'LEOP',              +0, 'LEOPActiveX',...
             'met2',        'dispEngine',   'LEOP',              +0, 'MET2',...
             'elsym5',      'dispEngine',   'LEOP',              +0, 'ELSYM5',...
             'nonlinear',   'bLinear',      1,                   +0, 0,...
             'igeophone',   'idxGeoPhone',  defIdxGeoPhone,      +1, [],...
             'resultasseed','bPrevAsSeed',  0,                   +0, 1,...
             'getpara',     'bGetPara',     0,                   +0, 1,...
             'noplot',      'bPlot',        1,                   +0, 0,...
             'help',        'NeedHelp',     0,                   +0, 1,...
             };

% Parse all of the options and saved into Options
[Options, sSpecified, unused] = ParseOptions(optionSet,varargin{:});
if Options.NeedHelp
    % help has already been shown
    return;
end

for i=1:length(unused)
    fprintf('WARNING: Unparsed input parameters - %s\n',num2str(unused{i}));
end

if any(Options.E==0)
    Message(0,1,'Dude, 0 is not a good guuess for layer stiffness! Try again man.');
    return
end

h = Options.h(:);
E = Options.E(:);
guessNMET2 = Options.guessNMET2;
if Options.bLinear
    XMap.bHasNMET2 = 0;
else
    XMap.bHasNMET2 = 1;
    if ~any(strcmp(Options.dispEngine,{'MET2','LEOP'}))
        Message(0,0,'KalmanBack: Only MET2 or LEOP can do non-linear back calculation\n');
        return;
    end
end

iDrop = Options.iDrop;
nStep = Options.nStep;
iIgnore = Options.iIgnore;
bResultAsSeed = Options.bPrevAsSeed;
idxGeoPhone = Options.idxGeoPhone;
iFixedE = Options.iFixedE;
iVarH = Options.iVarH(:);
dispEngine = Options.dispEngine;
bGetPara = Options.bGetPara;
XYLoad = Options.XYLoad;

% Start the ActiveX server, this will make it run way faster
% LEOP = actxserver('LEOPActiveX.LEOP');

if ~isempty(Options.eRange)
    XMap.eRange = UC('MPa','psi',Options.eRange);
end

% Allow some drops to be skipped
if isfield(Info,'iGroup')
    iGroup = Info.iGroup;
    iDrop = 1:size(iGroup,1);
else
    iGroup = 1:length(Drops);
    iGroup = iGroup(:);
end

iDrop = setdiff(iDrop,iIgnore);

% Determine XMap based on inputs
nLayer = length(h);    
XMap.idxE = setdiff(1:nLayer,iFixedE);
XMap.idxH = intersect(1:nLayer,iVarH);
if ~isempty(XMap.idxH)
    hRange = h(XMap.idxH)*[1,0.5];
    XMap.hRange = UC('mm','in',hRange);
end

switch Options.dispEngine
    case 'MET2'
        aDispServer = actxserver('LEOPActiveX.MET2');
    case 'LEOPActiveX'
        aDispServer = actxserver('LEOPActiveX.LEOP');
    otherwise
        aDispServer = '';
end
Options.DispServer = aDispServer;
    
m = 0;
for i=iDrop
    
    iOneGroup = iGroup(i,:);
    iOneGroup = iOneGroup(iOneGroup>0);
    
    sStatus = sprintf('Doing %5d/%5d',i,length(Drops));
    fprintf('%s\n',sStatus);
    
    para = [];
    disp = [];
    YMap = {};
    for j=iOneGroup
        temp = GetParam4EBack(Info,Drops(j),Options,unused);
        para = [para,temp.Para];
        disp = [disp;temp.Disp]; 
        YMap{end+1} = temp.YMap;
    end
    
    if bGetPara
        A = para;
        E0 = 0;
        return;
    end
    
    try 
        B = EKF_Main('para',para,'disp',disp,'engine',Options.dispEngine, ...
                 'guessE',E(XMap.idxE),'guessN',guessNMET2,'guessH',h(XMap.idxH), ...
                 'XMap',XMap,'YMap',YMap, ...
                 'TOL',0.05,unused{:},'nStep',nStep,'bPlot',Options.bPlot);    
    catch anError
        getReport(anError);
        Message(0,1,'EKF_Main: something wrong for drop %d\n',i);  
        continue;
    end
    
    % Post Processing the results
    B = MergeStruct(Drops(i),B);
    B.Ehist = UC('psi','mpa',B.State.E);
    if XMap.bHasNMET2
        B.nhist = B.State.nMET2;
        B.n = B.nhist(end);
    end    

    B.StationNo = Drops(i).PointNo;
    B.Station = Drops(i).Station;
    for j=1:6
        sEname = ['E',num2str(j)];
        if j<=nLayer
            val = B.Ehist(j,end);
        else
            val = B.Ehist(nLayer,end);
        end
        
        B.(sEname) = val;
    end
    
    hHist = UC('in','mm',B.State.h);
    for j=1:6
        sEname = ['h',num2str(j)];
        if j<=nLayer
            val = hHist(j,end);
        else
            val = hHist(nLayer,end);
        end
        
        B.(sEname) = val;
    end

    B.File = Drops(i).Filename;
    B.DropID = Drops(i).DropID;
    B.Time = Drops(i).Time;
    
    if isfield(Drops(i),'AirTemperature')
        B.AirTemperature = Drops(i).AirTemperature;
    end
    if isfield(Drops(i),'SurfaceTemperature')
        B.SurfaceTemperature = Drops(i).SurfaceTemperature;
    end
    if isfield(Drops(i),'Tpave')
        B.Tpave = Drops(i).Tpave;
    end
    if isfield(Drops(i),'Suspect')
        B.Suspect = Drops(i).Suspect;
    end
    
    B.Resample = Options.bResample;
    B.DispEngine = dispEngine;
    B.idxGeophone = idxGeoPhone;
    B.nGeophone = length(idxGeoPhone);
    B.Stress = Drops(i).Stress;
    B.rLoad = Drops(i).rLoad;
    B.ErrFinal = B.Err(end);
    B.nGuess = guessNMET2;
    B.guessE = E;
    
    if XMap.bHasNMET2
        para.E = B.State.E(:,end);
        para.nMET2 = B.n;
        [temp,B.C,B.Eom]=FunctMET2(para);
    end
    
    nGeoPhone = length(idxGeoPhone);
    B.RMS = sqrt(sum((B.Zhat-B.Z).^2)/nGeoPhone);
    pctErr = (B.Zhat./B.Z-1)*100;
    B.RMSPercent = sqrt(sum(pctErr.^2)/nGeoPhone);
    B.iDrop = i;
    
    B.Date = Info.Date;

    m = m+1;
    A(m) = B;
    
    
    xGeoPhone = UC('in','mm',para(1).XYout(:,1));
    xGeoPhone(xGeoPhone==0) = 75;
    
    nZout = para(1).nZout;
    nSamples = length(xGeoPhone);
    bSingle = nZout==1 & length(B.Z)==nSamples;
    if bSingle
        subplot(2,1,1)
    end
    
    plot(xGeoPhone,reshape(B.Z,[nSamples,length(B.Z)/nSamples]),'o'); hold on;
    plot(xGeoPhone,reshape(B.Zhat(:,end),[nSamples,size(B.Zhat,1)/nSamples]),'-'); hold off;
    
    if bSingle 
        pLoad = Drops(iOneGroup).Stress/1000;
        rLoad = Drops(iOneGroup).rLoad;
        Eot = (1-0.35^2)*pLoad*rLoad*rLoad./xGeoPhone;
        subplot(2,1,2);
        plot(Eot./B.Z,-xGeoPhone,'o',Eot./B.Zhat,-xGeoPhone,'r-');
        Eo.x = Eot./B.Z;
        Eo.y = xGeoPhone;        
    end
    title(sStatus);
    pause(0.1);
    
    if bResultAsSeed
        E = B.Ehist(:,end);
    end
end

if ~isempty(Options.DispServer)
    delete(Options.DispServer);
end

end


function [A] = GetParam4EBack(Info, Drop, Options, unused)

% stress must be in kPa
pLoad = Drop.Stress/1000;
rLoad = Drop.rLoad;

% the depth deflections are calculated
Zout = Drop.Depth;
disp = Drop.Def(:)/1000;        

if Options.bResample
    X = Info.rGeoPhone;
    Y = disp;    
    X1 = [0,200,300:Options.dSample:1500];
    Info.rGeoPhone = X1;
    disp = interp1(X,Y,X1,'spline');
    disp = disp(:);
    nGeoPhone = length(Info.rGeoPhone);
    Options.idxGeoPhone = 1:nGeoPhone;
end
%
% Determine which geophones to use
idxGeoPhone = Options.idxGeoPhone;
geoDist = Info.rGeoPhone(idxGeoPhone);
nGeoPhone = length(idxGeoPhone);
XYout = geoDist(:)*[1,0];

% Map the calculated deflections to the measured ones based on the deflection type
% This accounts for the reference points in MDD

switch upper(Drop.sDefType)
case {'FWD','RSD','EXX','EYY','EYY.M'},
    nZout = 1;
    nDisp = length(disp);
    nTtlGeoPhone = length(Info.rGeoPhone);        
    nGroup = nDisp/nTtlGeoPhone;
    YMap = eye(nGeoPhone);
    YMap = repmat(YMap,[nGroup,1]);
    
    idx = zeros(nTtlGeoPhone,1);
    idx(idxGeoPhone) = 1;
    idx = repmat(idx,[nGroup,1]);        
    disp = disp(idx>0);
    
    XYLoad = [0,0];
    if upper(Drop.sDefType(1))=='E'
        XYout(:,2) = Drop.yOut;
        XYLoad = Info.XYLoad;
        Zout = Drop.Depth;
        disp = UC('in','mm',disp)*1000; % to cancel the transformation
    end           
    
case 'MDD',
    Zout = [Zout,3000];
    nZout = length(Zout);
    nDisp = nGeoPhone*nZout;
    YMap = eye(nDisp-nGeoPhone);
    YMap = [YMap, repmat(-eye(nGeoPhone),[nZout-1,1])];
case {'RSD+MDD','MDD+RSD'},
    Zout = [Zout,3000];
    nZout = length(Zout);
    nDisp = nGeoPhone*nZout;
    YMap = eye(nDisp-nGeoPhone);
    
    YMap = [YMap, repmat(-eye(nGeoPhone),[nZout-1,1])];
    YMap(1:nGeoPhone,nDisp-nGeoPhone+1:end) = 0;                   
    
otherwise,
    fprintf('Unknown Deflection Type: %s\n',Info.sDefType);
    return;
end

para = ELSYM5_Init(Options.dispEngine,...
                   'layer',Options.h(:),'v',0.35,'E',Options.E(:), ...
                   'pressure',pLoad,'rLoad',rLoad,'xyLoad',XYLoad, ...                       
                   Drop.sDefType, ...
                   'XYout', XYout,'Zout',Zout, unused{:});

% transfer the DispServer field
para.DispServer = Options.DispServer;

A.Para = para;
A.Disp = disp;
A.YMap = YMap;
               
end     
