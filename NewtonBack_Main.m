function [A] = NewtonBack_Main

para = Elsym5_Init('MET2','e',[5000;400;100],'h',[150;400;0]);

para.XMap = struct('idxE',1:3, 'bHasNMET2', 0);
para.Z = DispEngine(para); % takes untransformed E, n

bounds = log(uc('MPa','psi',repmat([10 30000],[3,1])));
[topGenes, fitness] = Newton(@CostFunctDispLET, para);

A.topGenes = uc('psi','Mpa',exp(topGenes));
A.fitness = fitness;

