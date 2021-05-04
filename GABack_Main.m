function [A] = GABack_Main(para,XMap,zMeasure,bounds,nGenes,nIter)

%para = Elsym5_Init('MET2','e',[400;4000;100],'h',[150;400;0]);
%para.XMap = struct('idxE',1:3, 'bHasNMET2', 0);
%para.Z = DispEngine(para); % takes untransformed E, n
%bounds = log(uc('MPa','psi',repmat([10 30000],[3,1])));

[topGenes, fitness] = GeneticOptimize(bounds, @CostFunctDispLET, nGenes, nIter, para);

A.topGenes = uc('psi','Mpa',exp(topGenes));
A.fitness = fitness;

