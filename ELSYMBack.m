function [Eback,para,fitness]=ELSYMBack(E,data)

layers = data.Layers;
n = length(layers);
% Combine the AB/ASB layers
%    layers = [layers(1:n-3);layers(n-2)+layers(n-1);layers(n)];
if(n==5)
    layers = [sum(layers(1:2));sum(layers(3:4));layers(end)];
else
    layers = [sum(layers(1:3));sum(layers(4:5));layers(end)];
end

nLayers = length(layers);
para = struct('Load',Load,'v',0.35*ones(1,nLayers));

switch type
case 'MDD',
    MDD = data;
    
    data = MDD.Data;
    flt = MyStruct('ID',Station,'Rep',Rep,'Load',Load);
    anchor = MDD.Anchor;
    
    idx = ApplyFilter(data,flt);
    def = [[data(idx).Depth];[data(idx).Def]]';
    def = sortrows(def,1);
    
    dDef = diff(def(:,2));
    if(any(dDef>0))
        Eback = -1;
        return;   
    end

    def = [def(:,1)*0,def];
    
    rOffset = UnitConvert('mm','in',200);  % for dual wheel tire, set to 0 for single tire
    rMea = unique(def(:,1));
    zMea = [unique(def(:,2));anchor];
    zDef = def(:,3);
    
    pTire = 0.69;
    
case 'FWD',
    FWD = data;
    rMea = FWD.rMea;
    zMea = 0;
    rOffset = 0;
    zDef = mean(FWD.da)'/1000;
    rLoad = 150;
    pTire = Load/2*1000/pi/rLoad/rLoad;
end
    
    
nLayers = length(layers);

% Make sure the number of layers are within limits

if(nLayers>5)
    fprintf('ERROR: No more than 5 layers (include subgrade) allowed\n');
    return;
end

% Convert the units into inch, if necessary
if(max(layers)<30)
    fprintf('FYI: Assume you are using inches for layer thickness\n');
    bInchIn = 1;
else
    fprintf('FYI: Assume you are using mm for layer thickness\n');
    layers = UnitConvert('mm','in',layers);
    rMea = UnitConvert('mm','in',rMea);
    zMea = UnitConvert('mm','in',zMea);
    zDef = UnitConvert('mm','in',zDef);
    bInchIn = 0;
end
layers(end) = 20*12-sum(layers);

if(pTire<10)
    fprintf('FYI: Assume you are using MPa for tire pressure\n');
    pTire = UnitConvert('MPa','psi',pTire);
else
    fprintf('FYI: Assume you are using psi for tire pressure\n');
end    

if(Load<500)
    fprintf('FYI: Assume you are using kN for Axle load\n');
    Load = UnitConvert('kN','lbs',Load);
else
    fprintf('FYI: Assume you are using lbs for Axle load\n');
end

para.zDef = zDef;

% Setup the necessary input for ELSYM5 and Generate the input file
para.title='Back Calculation for layer stiffnesses';
para.Load=Load/2;
para.Pressure=pTire;
radius = sqrt(para.Load/pi/para.Pressure);
para.nXYout=length(rMea);
para.XYout=(rMea(:)+rOffset)*[1,0];
para.nLayers=length(layers);
%para.E=para.E;
para.h=layers;
para.v=para.v;
para.nZout=length(zMea);
para.Zout=zMea;

% For cases without ATPB, which assume 1000MPa for the stiffness
if(nLayers==4)
    Ebounds = [1000,10000;100,5000;100,1000;10,500] %;10,500]; % in MPa;
elseif(nLayers==5)
    % For caes with ATPB at the 3rd layer
    Ebounds = [1000,10000;100,5000;1000,1000;100,1000;10,500] %;10,500]; % in MPa;
elseif(nLayers==3)
    Ebounds = [100,10000;100,1000;10,500]
end

Ebounds = log(UnitConvert('MPa','psi',Ebounds));

[Eback, fitness] = GeneticOptimize(Ebounds,@FunctElsym5,100,21,para);

Eback = UnitConvert('psi','MPa',exp(Eback));




