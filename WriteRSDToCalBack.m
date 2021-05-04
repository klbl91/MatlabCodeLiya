function [E,A] = WriteRSDToCalBack(file,data,sSection)

nRSD = length(data);
loadLevels = unique([data.TestLoad]);

% Copy the target database as the dummy DB
PrepareODBC(file);

for i=1:length(loadLevels)
    
    idx = find([data.TestLoad]==loadLevels(i));
    Reps(i).Value = unique([data(idx).Repetition]);
    
    m = min(idx);
    RSD = data(m);
    
    fprintf('Doing Load level %d\n',RSD.TestLoad);

    Force = RSD.TestLoad/2; % kN;
    p = 0.72; % MPa;
    rLoad = sqrt(Force*1000/p/pi); % mm
    
    % using equivalent single wheel load
    geoPhone = [0 200 300 460 610 910 1525];
    geoPhone = sqrt(geoPhone.^2+150^2); % Assuming the two wheels are 300-mm apart

    % First Table "Geo"
    
    % Title is identifier for a set of data with the same geophone settings
    sTitle = [sSection,'_',num2str(RSD.TestLoad),'kN'];
    
    TableName = 'Geo';
    ColNames = {'Filename','Numofgeophones','RadiusOfPlate'};
    NewData = {sTitle, 7, rLoad};
    for j=1:7
        ColNames{end+1} = ['DistFrCntr',num2str(j)];
        NewData{end+1} = geoPhone(j);
    end
    
    WriteDatabase(GetDummyDB('DSN'),TableName,ColNames,NewData);

    % Then Table "Fwdfiles"
    TableName = 'Fwdfiles';
    ColNames = {'Filename','"Date"','"Road name"','Direction','Lane','"From"','"to"','"Road Id"'};
    NewData = {sTitle,'20060601','x','x',0, 1, length(Reps(i).Value), 0};
    
    WriteDatabase(GetDummyDB('DSN'),TableName,ColNames,NewData);

end


% Deflections
TableName = 'Deflections';
ColNames = {'Filename','Chainage','Temperature1','"Drop No"','Stress', ...
            'D1','D2','D3','D4','D5','D6','D7','"Test point"','Position','TheTime','GPS','"Note"'};
    
clear NewData;

NewDate = {sSection};
iDropMap = [[2000,3000,4000,5000,6000]+500,[2000,4000,6000]+300,[2000,4000,6000]+700];

for i=1:nRSD
    RSD = data(i);
    
    Force = RSD.TestLoad/2; % kN;
    p = 0.72; % MPa;
    rLoad = sqrt(Force*1000/p/pi); % mm
    
    NewData{i,1} = [sSection,'_',num2str(RSD.TestLoad),'kN'];
    NewData{i,2} = RSD.Repetition;
    NewData{i,3} = UnitConvert('c','f',20);
    
    if RSD.Y<500
        RSD.Y = 300;
    elseif RSD.Y>500
        RSD.Y = 700;
    end
    iDrop = find(iDropMap==RSD.X+RSD.Y);                
    NewData{i,4} = iDrop;
    
    NewData{i,5} = p*2*1000; % in kPa
    for j=1:7
        NewData{i,j+5} = RSD.yGeoPhone(j)*1000; % in micron
    end
    
    iLoad = find(loadLevels==RSD.TestLoad);
    m = find(Reps(iLoad).Value==RSD.Repetition);
    NewData{i,13} = m;
    NewData{i,14} = 1;
    NewData{i,15} = 1001;
    NewData{i,16} = 'NA';
    NewData{i,17} = NewData{i,1};
    
end

iDrops = [NewData{:,4}];
idx = find(iDrops<=5);
NewData = NewData(idx,:);
WriteDatabase(GetDummyDB('DSN'),TableName,ColNames,NewData);

% Sectioning

