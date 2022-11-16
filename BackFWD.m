function [A,Eo] = BackFWD(A,varargin)

fprintf('%s\n',A.File);

Drops = A.Drops;
Info = A.Info;

clear A;

% Default parameters

if isfield(Info,'Thickness')
    layers = Info.Thickness(:);
else
    layers = [125; 349; 0];
end

guessE = [1000;200;80];
guessNMET2 = -0.2;
XMap.bHasNMET2 = 1;
iDrop = 1:length(Drops);
iIgnore = [];
bResultAsSeed = 0;
idxGeoPhone = [];
iFixed = [];
dispEngine = 'MET2';
bGetPara = 0;

i = 1;
j = 0;
args = {};
while i<=length(varargin)
    sOption = varargin{i};
    switch lower(sOption)
    case 'n',
        guessNMET2 = varargin{i+1};
        i = i+1;
    case 'idrop',
        iDrop = varargin{i+1};
        i = i+1;
    case 'ignore',
        iIgnore = varargin{i+1};
        i = i+1;
    case 'h',
        layers = varargin{i+1};
        layers = layers(:);
        i = i+1;
    case 'e',
        guessE = varargin{i+1};
        guessE = guessE(:);
        i = i+1;
    case 'fixed',
        iFixed = varargin{i+1};
        iFixed = iFixed(:)';
        i = i+1;
    case {'linear','linear.met2'},        
        XMap.bHasNMET2 = 0;
        dispEngine = 'MET2';
    case {'linear.elsym5','elsym5'},
        XMap.bHasNMET2 = 0;
        dispEngine = 'ELSYM5';
    case 'nonlinear',
        XMap.bHasNMET2 = 1;
        dispEngine = 'MET2';        
    case 'igeophone',
        idxGeoPhone = varargin{i+1};
        i = i+1;
    case 'resultasseed',
        bResultAsSeed = 1;
    case 'getpara',
        bGetPara = 1;
    otherwise,
        j = j+1;
        args{j} = varargin{i};
    end
    i = i+1;
end


if length(idxGeoPhone)==0
    idxGeoPhone = 1:length(Info.rGeoPhone);
end
geoDist = Info.rGeoPhone(idxGeoPhone);
XYout = geoDist(:)*[1,0];

m = 0;
iDrop = setdiff(iDrop,iIgnore);

for i=iDrop
    
    fprintf('Doing %5d/%5d\n',i,length(Drops));
    pLoad = Drops(i).Stress/1000;
    if isfield(Drops(i),'rLoad')
        rLoad = Drops(i).rLoad;
    else
        rLoad = 150; %sqrt(F/pi/pLoad);
    end
    para = Elsym5_Init(dispEngine, Info.sDefType, ...
                       'layer',layers,'v',0.35,'E',guessE, ...
                       'pressure',pLoad,'rLoad',rLoad, ...
                       'XYout', XYout,args{:});
    
    
    if bGetPara
        A = para;
        E0 = 0;
        return;
    end
    
    k = 0;
    for j=idxGeoPhone
        k = k+1;
        disp(k) = getfield(Drops(i),['D',num2str(j)]);
    end
    disp = disp(:)/1000;
    
    nLayer = length(layers);
    
    XMap.idxE = 1:nLayer;
    XMap.idxE = setdiff(XMap.idxE,iFixed);
    guessX = [log(uc('mpa','psi',guessE(XMap.idxE)));tan(guessNMET2*pi/2)];
    
    
    B = EKF_Main('para',para,'disp',disp,'guess',guessX,'XMap',XMap,'TOL',0.1,args{:});    
    
    
    B.Ehist = guessE*ones(1,size(B.X,2));
    B.Ehist(XMap.idxE,:) = uc('psi','mpa',exp(B.X(1:end-1,:)));
    if XMap.bHasNMET2
        B.nhist = atan(B.X(end,:))/(pi/2);
        B.n = B.nhist(end);
    end    

    for j=1:nLayer
        sEname = ['E',num2str(j)];
        B = setfield(B,sEname,B.Ehist(j,end));
    end
    B.Station = Drops(i).Station;
    B.DropID = Drops(i).DropID;
    if isfield(Drops(i),'AirTemperature')
        B.AirTemperature = Drops(i).AirTemperature;
    end
    if isfield(Drops(i),'SurfaceTemperature')
        B.SurfaceTemperature = Drops(i).SurfaceTemperature;
    end
    B.Stress = Drops(i).Stress;
    B.ErrFinal = B.Err(end);
    B.nGuess = guessNMET2;
    B.guessE = guessE;
    B.h = layers;
    
    if XMap.bHasNMET2
        para.E = uc('MPa','psi',[B.E1;B.E2;B.E3]);
        para.nMET2 = B.n;
        [temp,B.C,B.Eom]=FunctMET2(para);
    end
    
    nGeoPhone = length(idxGeoPhone);
    B.RMS = sqrt(sum((B.Zhat-B.Z).^2)/nGeoPhone);
    B.RMSPercent = sqrt(sum((B.Zhat./B.Z-1).^2)/nGeoPhone);
    
    E(:,i) = B.Ehist(:,end);
    
    B.iDrop = i;
    
    m = m+1;
    A(m) = B;
    
    
    xGeoPhone = uc('in','mm',para.XYout(:,1));
    xGeoPhone(xGeoPhone==0) = 75;
    
    subplot(2,1,1)
    plot(xGeoPhone,B.Z,'o',xGeoPhone,B.Zhat,'r-');
    
    Eot = (1-0.35^2)*pLoad*rLoad*rLoad./xGeoPhone;
    subplot(2,1,2);
    plot(Eot./B.Z,-xGeoPhone,'o',Eot./B.Zhat,-xGeoPhone,'r-');
    Eo.x = Eot./B.Z;
    Eo.y = xGeoPhone;
    pause(0.1);
    
    if bResultAsSeed
        guessE = B.Ehist(:,end);
    end
end
