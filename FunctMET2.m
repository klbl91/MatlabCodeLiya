function [A,C,Eom]=FunctMET2(para)

nProb = size(para.E,2);

F = para.Load;
pLoad = para.Pressure;
rLoad = sqrt(F/pLoad/pi);
nGeophone = para.nXYout;

aServer = para.DispServer;
if ~isempty(aServer)
    MET2 = aServer;
    bServerCreated = false;
else
    MET2 = actxserver('LEOPActiveX.MET2');
    bServerCreated = true;
end

set(MET2,'Parameters','rLoad',UC('in','mm',rLoad));
set(MET2,'Parameters','pLoad',UC('psi','mpa',pLoad));

set(MET2,'Parameters','nLayer',para.nLayers);
set(MET2,'Parameters','nGeoPhone',nGeophone);
set(MET2,'Parameters','Rgeophone',UC('in','mm',para.XYout(:,1)));

C(nProb) = 0;
disp(1:nGeophone,nProb) = 0;
Eom(nProb) = 0;

for i=1:nProb
    set(MET2,'Parameters','Stiffness',UC('psi','mpa',para.E(:,i)));
    set(MET2,'Parameters','Thickness',UC('in','mm',para.h(:,i)));
    set(MET2,'Parameters','n',para.nMET2(i));
    invoke(MET2,'CallMET2');
    temp = get(MET2,'Deflections');
    disp(:,i) = temp(:);
    C(i) = get(MET2,'Results','c');
    Eom(i) = get(MET2,'Results','Eom');
end

if bServerCreated
    delete(MET2);
end

A = UC('mm','in',disp);
