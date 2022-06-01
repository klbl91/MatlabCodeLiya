%%4/12/2021
%%LJ
%%------------SCB parameters calculation------------%
%--- codes to calculate the SCB paramters from load-displ curves of brittle fracture----%
% Input: SCB geometries (Ligament, Thickness, Specimen ID, Radius);
%       SCB testing results (load-displ curves);
% Ouput: SCB parameters (Fracture engery, FI, FIasc, front slope,...
%       back slope, KIC, strength, prepeak area, afterpeak area)
% Adopted from CalculatedSCBParameters_LJ.m
% The load vs displacement of brittle fracture couldn't be fitted with Gauss4 due to 
%   limitation of testing points (9 points minimum required for Gauss4 fit)
%   Gauss2->1 is used here for fitting and the inflection point is the local minimum 
%   of second derivative of Gauss1 fitting result 
%------------------------------------------------------------------
% 5/31/2022
% modification: output the .DB file for OLTS

%% calculation steps need to update for Index=1 and Index=0 following

clc;
clear;
testdatapath=pwd;

[basename,folder]=uigetfile(testdatapath+"\*.xlsx",'select the input data file');
fullFileName=fullfile(folder,basename);

s=importdata(fullFileName);

A='Ligament';%i 
B='Thickness';%j
C='Specimen ID';%k
D='Radius';%r
for i=1:size(s.textdata,2)
    if contains(s.textdata{1,i},A)
        break       
    end
    if i==size(s.textdata,2)&&~contains(s.textdata{1,i},A)
       msg='Error: Could not locate the information: '+convertCharsToStrings(A);
       error(msg);
    end
end
for j=1:size(s.textdata,2)
    if contains(s.textdata{1,j},B)
        break
    end
    if j==size(s.textdata,2)&&~contains(s.textdata{1,j},B)
       msg='Error: Could not locate the information: '+convertCharsToStrings(B);
       error(msg);
    end
end
for k=1:size(s.textdata,2)
    if contains(s.textdata{1,k},C)
        break
    end
    if k==size(s.textdata,2)&&~contains(s.textdata{1,k},C)
       msg='Error: Could not locate the information: '+convertCharsToStrings(C);
       error(msg);
    end
end
for r=1:size(s.textdata,2)
    if contains(s.textdata{1,r},D)
        break
    end
    if r==size(s.textdata,2)&&~contains(s.textdata{1,r},D)
       msg='Error: Could not locate the information: '+convertCharsToStrings(D);
       r=75;
    end
end
disp("The specimens that need to be calculated are:")
disp(s.textdata(2:end,k));

[basename2,folder2]=uigetfile(testdatapath+"\*.csv",'select the SCB testing results file (s)','Multiselect','on');
fullFileName2=fullfile(folder2,basename2);
if ~iscell(fullFileName2)
    fullFileName2={fullFileName2};
    basename2={basename2};
end

if size(fullFileName2,2)~=size(s.data,1)
    msg="The testing files number does not match with the input file samples number";
    error(msg);
end

sFields = {
    'Name',...
    'Testduration_',...
    'MaxLoadkN','Strengthpsi','ResidualLoadkN',...
    'SecantModuluskNmm',...    
    'DisplacementChannel1Rate_',...
    'MaxDisplacementChannel1mm_',...
    'DisplacementChannel1atFailuremm_',...
    'DisplacementChannel1basedFract1',...
    'DisplacementChannel1basedFract2',...
    'PrepeakareaDisplacementChannel1',...
    'InterceptDisplacementChannel1_',...
    'SlopeDisplacementChannel1_',...
    'FlexibilityIndexDisplacementCh',...
    'InitialDisplacementTrim_',...
    'ResidualLoadTrim_',...
    'SpecimenID',...
    'Date',...
    'DisplacementRate',...
    'FlexibilityIndex',...
    'Comments',...
    'Gf','Sigma_p','Spp','TestMachine'
    %'#Inflection point found at',... %    disp 1.667180e+00 mm and load 3.853632e+00 kN.
};
sFields = sFields';
nFields = length(sFields);


for m=1:size(fullFileName2,2)
    DATA=xlsread(fullFileName2{m});
    TestsampleID=extractBefore(basename2{m},' ');
    Testdate=extractBetween(basename2{m},' ','.csv');
    Datetime=datetime(Testdate{1},'InputFormat','MM-dd-yy hh.mm.ss a',...
        'Format','MM/dd/uuuu hh:mm:ss aa');
    Testduration=max(DATA(:,2));%sec
    [peakload,peakrow]=max(DATA(:,3));%max load (kN)
    Strength=peakload*10^9/(2*s.data(m,r)*s.data(m,j));%%unit N/m2
    StrengthPSI=Strength*0.000145;
    ResidualLoad=DATA(end,3);%kN
    preU=DATA(1:peakrow,4);
    preLoad=DATA(1:peakrow,3);
    postU=DATA(peakrow:length(DATA),4);
    postLoad=DATA(peakrow:length(DATA),3);
    frontarea=trapz(preU,preLoad);
    backarea=trapz(postU,postLoad);
    DisplacementRate=(DATA(end,4)-DATA(2,4))/(DATA(end,2)-DATA(2,2));
    MaxDisplacement=DATA(end,4);
    DisplacementatFailure=DATA(peakrow,3);
    DisplacementbasedFracturearea_short=inf;%not calculated
    DisplacementbasedFracturearea_long=inf;%not calculated
    Intercept=inf;

    %%calculate Gf
    %fitting of curve before peak
    p1=polyfit(preU,preLoad,6);
    q1=polyint(p1);
    Wf1=diff(polyval(q1,[0,DATA(peakrow,4)]));
    preloadclean=polyval(p1,preU);
    %fitting of curve after peak
    try 
        p2=fit(postU,postLoad,'gauss2');
    catch 
        warning("fitting issue with gauss2, continue to use gauss2")
        p2=fit(postU,postLoad,'gauss1');
    end
    postloadclean=p2(postU);

    
    %DispEnd=arrayfun(@(y)fzero(@(x)p2(x)-y,1),0.1);%calculate the corresponding displacement when loading =0.1kN
    [~,INDEXdispEnd]=min(abs(DATA(peakrow:length(DATA),3)-0.1));
    DispEnd=DATA(INDEXdispEnd+peakrow-1,4);
    Wf2=abs(integrate(p2,DATA(peakrow,4),DispEnd));
    Wf=Wf1+Wf2;
    Gf=Wf*10^6/((s.data(m,i))*s.data(m,j));
    %%calculate slope
    DATA1=DATA(1:peakrow,:);
    [~,point1r]=min(abs(DATA1(:,3)-1/4*DATA1(peakrow,3)));
    [~,point2r]=min(abs(DATA1(:,3)-3/4*DATA1(peakrow,3)));
    Sasc=(DATA1(point2r,3)-DATA1(point1r,3))/(DATA1(point2r,4)-DATA1(point1r,4));%kN/mm
    Secantmodulus=Sasc;

    PeakDisp=DATA(peakrow,4);
    tt=1;
    for incr=PeakDisp:0.0001:DispEnd
        [~,inflect(tt)]=(differentiate(p2,incr));
        tt=tt+1;
    end
    localMin=islocalmin(abs(inflect));
    inflectIndex=find(localMin==1);
    inflectIndex=inflectIndex(1);
    inflect_pt=PeakDisp+0.0001*inflectIndex;
    [~,inflect_ptindex]=min(abs(postU-inflect_pt));
    inflect_load=postLoad(inflect_ptindex);
    Spp=double(differentiate(p2,inflect_pt));
    %%calculate the FI
    FI=Gf*0.01/abs(Spp);
    %FIasc=Gf*0.01/Sasc;
    %%calculate KIC
    %YI=4.782+1.219*(s.data(m,r)-s.data(m,i))/s.data(m,r)+0.063*exp(7.045*(s.data(m,r)-s.data(m,i))/s.data(m,r));    %TT.YI=4.782+1.219*(TT.Notch./TT.Radius)+0.063*exp(7.045*(TT.Notch./TT.Radius));
%TT.KICmatlab=TT.YI.*sqrt(pi()*TT.Notch/1000).*(TT.Pmax*1e3./(2*TT.Radius.*TT.TT.Thickness);Notch,Radius and thickness are in mm; Pmax in kN                                                                    

    %KIC{m}=YI*sqrt(pi()*(s.data(m,r)-s.data(m,i))/1000)*peakload*10^3/(2*s.data(m,r)*s.data(m,j));%the unit of KIC is MPa(sqrt(m))
            
       
    figure
    plot(DATA(:,4),DATA(:,3));
    title(TestsampleID);
    xlabel('Deformation (mm)');
    ylabel('Load (kN)');
    hold on
    syms xx
    y = Spp*(xx-inflect_pt)+inflect_load;
    fplot(y, [PeakDisp DispEnd], '--r')
    plot(inflect_pt,inflect_load,'r*')
    ylim([0 peakload])
    hold off
    saveas(gcf,strcat(folder,TestsampleID,'.png'))
    
 %%Fracture engery, FI, FIasc, front slope,
%%       back slope, KIC, strength, prepeak area, afterpeak area)

Outputsummary={'Value',...
    Testduration,...
    peakload,...
    StrengthPSI,...
    ResidualLoad,Secantmodulus,DisplacementRate,MaxDisplacement,DisplacementatFailure,...
    "NONE","NONE",Wf1,"NONE",Spp,FI,...
    "No trimming",0.1,...
   TestsampleID,Datetime, DisplacementRate,FI,...
    "High brittleness",Gf,StrengthPSI,Spp,"SCB#1"
    };
loadclean=[preloadclean(1:end-1);postloadclean];
cleandata=[DATA(:,2),loadclean,DATA(:,4)];
cleandata=array2table(cleandata);
cleandata.Properties.VariableNames={'TimeSec','LoadkN','Displacementmm'};

Outputsummary=Outputsummary';
writecell(sFields,strcat(folder,TestsampleID,"_brittle_DB.xlsx"),'Sheet','Data.Summary','Range','A1');
writecell(Outputsummary,strcat(folder,TestsampleID,"_brittle_DB.xlsx"),'Sheet','Data.Summary','Range','B1');

writetable(cleandata,strcat(folder,TestsampleID,"_brittle_DB.xlsx"),'Sheet','Date.Clean','WriteVariableNames',true);
end
