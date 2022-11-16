function [A]=TestMETActiveX

pLoad = 0.934;
rLoad = 150;
geoPhone = [0 200 300 600 900 1200 1500];
nLayers = 3;
h = [100;400;0];
E = [5000;300;100];

% initialize the ActiveX server, you need to run MET2ActiveX.exe at least once
% to register this server

%MET2 = actxserver('MET2ActiveX.MET2');
WorMET2 = actxserver('LEOPActiveX.MET2');

% set up the parameters, the parameter name is NOT case-sensitive
% all the units are in MPa and mm

set(MET2,'Parameters','rLoad',rLoad);
set(MET2,'Parameters','pLoad',pLoad);
set(MET2,'Parameters','nLayer',nLayers);
set(MET2,'Parameters','Thickness',h);
set(MET2,'Parameters','Stiffness',E);
set(MET2,'Parameters','nGeoPhone',length(geoPhone));
get(MET2,'Parameters','nGeoPhone');
set(MET2,'Parameters','Rgeophone',geoPhone);
set(MET2,'Parameters','n',-0.3);

% You can retrieve the parameter values like this
A.xGeoPhone = get(MET2,'Parameters','Rgeophone');

% now ask the activeX server to do the calculation
invoke(MET2,'CallMET2');

% then you can get the results
A.yGeoPhone = get(MET2,'Results','Deflections');
A.c = get(MET2,'Results','c');
A.Eom = get(MET2,'Results','Eom');

% after you are done, release the server
delete(MET2);

