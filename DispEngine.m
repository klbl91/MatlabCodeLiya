function [disp]=DispEngine(para)

switch para.DispEngine
    case 'ELSYM5'
        disp = FunctELSYM5(para);
    case 'MET2'
        disp = FunctMET2(para);
    case 'LEOP'
        disp = FunctLEOP(para);
    case 'LEOPActiveX'
        disp = FunctLEOPActiveX(para);
end

switch para.sDefType
case 'RSD'
    OA = UC('mm','in',500);
    OB = UC('mm','in',3000);
    AB = OB-OA;
    R0 = para.XYout(1,1);
    rGeoPhone = sqrt(para.XYout(:,1).^2-R0^2);
    
    paraTmp = para;
    paraTmp.sDefType = 'FWD';    
    
    paraTmp.XYout(:,1) = sqrt((rGeoPhone+OB).^2+R0^2);    
    dispO = DispEngine(paraTmp);
    
    paraTmp.XYout(:,1) = sqrt((rGeoPhone+AB).^2+R0^2);    
    dispA = DispEngine(paraTmp);
    
    K = OB/OA;
    % dMea = dB+(K-1)*dO-K*dA
    %disp = (disp-dispO)+K*(dispA-dispO);
    disp = disp+(K-1)*dispO-K*dispA;
case 'EYY.M'
    paraTmp = para;
    paraTmp.sDefType = 'EYY';

    % Note: need to adjust this 4000 later for the actual staring distance
    % to the strain gauge
    paraTmp.XYout(end+1,:) = [4000/25.4 0];
    paraTmp.nXYout = para.nXYout+1;
    disp = DispEngine(paraTmp);
    nDisp = size(disp,1)-1;
    disp = disp(1:nDisp,:)-repmat(disp(end,:),nDisp,1);
end

if any(para.bRelative)
    nDisp = size(disp,1)/para.nZout;
    dispRef = disp(end-nDisp+1:end,:);
    for i=1:(para.nZout-1)
        if ~para.bRelative(i)
            continue;
        end
        
        idx = (nDisp*(i-1)+1):nDisp*i;
        disp(idx,:) = disp(idx,:)-dispRef;
    end
    disp = disp(1:end-nDisp,:);
end
