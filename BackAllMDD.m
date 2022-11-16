function [E,A] = BackAllMDD(data)

nMDD = length(data);

for i=1:nMDD
    
    MDD = data(i);
    
    fprintf('Doing %10d%10d(%10d/%10d)\n',MDD.iMDD,MDD.Repetition,i,nMDD);

    Force = MDD.TestLoad/2; % kN;
    p = 0.72; % MPa;
    rLoad = sqrt(Force*1000/p/pi); % mm
    para = Elsym5_Init('layer',[125;349;0], ...
                        'v',0.35,'E',[8000;500;200], ...
                        'load',Force,'rLoad',rLoad, ...
                        'xyload',[0 150;0 -150], ...
                        'xyout',MDD.xGeoPhone(:)*[1 0], ...
                        'Zout',[MDD.zLVDT],'Relative',MDD.bRelative);

    disp = MDD.yGeoPhone;
    disp = disp(:);
    if length(disp)<3
        fprintf('No enough LVDT working\n');
        continue;
    end
    
    guessE = [2000;100;50];
    
    try
        A(i) = EKF_Main('para',para,'disp',disp,'guess',guessE);
        E(:,i) = A(i).X(:,end);
        bError(i) = 0;
        
    catch
        A(i).X = 0;
        E(1:3,i) = 0;
        fprintf('ERROR occured for i=%d\n',i);
        bError(i) = 1;
        fprintf('%s\n',lasterr);
    end

    save BackMDDResults E A bError
    
end

MDD = data;

save BackMDDResults MDD -APPEND;
