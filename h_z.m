function [Zhat]=h_z(X,Input_geometry,XMap,YMap)

Zhat = [];
for i=1:length(Input_geometry)
    % Map X onto structure parameters
    
    param = MapXtoState(X,XMap,Input_geometry(i));
    
    % Finally call the displacement engine
    
    disp = DispEngine(param);
    Zhat = [Zhat;YMap{i}*disp];
end