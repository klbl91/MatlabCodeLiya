function [A] = MapStateToX(E,nMET2,h,XMap)

if nargin==1
    Options = E;
    XMap = Options.XMap;
    E = Options.guessE;
    nMET2 = Options.guessNMET2;
    h = Options.guessH;
end

A = [];
if ~isempty(E)
    if isfield(XMap,'eRange')
        E = (E-XMap.eRange(:,1))./XMap.eRange(:,2);
        E = tan(E*pi/2);
    else
        E = log(E-10);
    end
    A = E;
end

if ~isempty(h)
    h = (h-XMap.hRange(:,1))./XMap.hRange(:,2);
    h = tan(h*pi/2);
    A = [A;h];
end

if XMap.bHasNMET2
    A = [A;tan(nMET2*pi/2)];
end

