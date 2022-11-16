function [A] = MapXtoState(X,XMap,Para)

nProb = max([size(X,2),1]);
nEinX = length(XMap.idxE);
nHinX = length(XMap.idxH);

if nEinX>0
    idxE = 1:nEinX;
    if isfield(XMap,'eRange')
        E = atan(X(idxE,:))/(pi/2);
        
        E = [E; ones(1,nProb)];
        scale = diag(XMap.eRange(:,2));
        offset = XMap.eRange(:,1);
        E = [scale,offset]*E;
    else
        E = exp(X(1:nEinX,:))+10;
    end
    State.E = E;
end

if nHinX>0
    idxH = nEinX+[1:nHinX];
    H = atan(X(idxH,:))/(pi/2);
    H = [H; ones(1,nProb)];
    hRange = XMap.hRange;
    scale = diag(hRange(:,2));
    offset = hRange(:,1);
    xTrans = [scale,offset];
    H = xTrans*H;
    State.H = H;
end
    
if XMap.bHasNMET2
    nMET2 = atan(X(end,:))/(pi/2);
    State.nMET2 = nMET2;
end

% Return the State vector in different ways

if nargin==3
    Para.E = Para.E*ones(1,nProb);
    Para.h = Para.h*ones(1,nProb);
    Para.nMET2 = Para.nMET2*ones(1,nProb);
    
    if nEinX>0
        Para.E(XMap.idxE,:) = E;
    end
    
    if nHinX>0
        Para.h(XMap.idxH,:) = H;
    end
    
    if XMap.bHasNMET2
        Para.nMET2 = nMET2;        
    end
    
    A = Para;

else    
    A = State;
end


