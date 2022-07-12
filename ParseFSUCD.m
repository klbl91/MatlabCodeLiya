function [] = ParseFSUCD()
%To run: ParseFSUCD
%Analyze one mixture a time
% Unless a frequency different from 10Hz:ParseFatUCD(0.05) FOR 0.05 Hz
% require the two output excel files from the UTS015
% require the empty excel processor "4PB Fatigue ProcessorforOLTS.xlsx"
% require function sinFitSSE
[basename,folder]=uigetfile(('*/*Runtime*.csv'),'select all the input data files for one mixture only(s)','Multiselect','on');%read raw data file
fullFileName=fullfile(folder,basename);
while length(unique(extractBefore(basename,"(FAT")))>1%check if it is a single mixture
 reselmsg=msgbox("You have selected multiple mixtures! Please reselect");
 [basename,folder]=uigetfile(('*/*Runtime*.csv'),'select all the input data files for one mixture only(s)','Multiselect','on');%read raw data file
 fullFileName=fullfile(folder,basename);
 delete (reselmsg);
end
MachineID=5;
if contains(strrep(lower(basename{1})," ",""),"fat6")
    MachineID=6;
end
% if ~iscell(fullFileName)
%     fullFileName={fullFileName};
%     basename={basename};
% end
%hw = waitbar(0,'Running...');
frequency={[]};
p2p_stress={[]};%MPa
p2p_strain={[]};%mm/mm
EG={[]};%MPa
PhaseAngle={[]};
Temperature={[]};
% specimen dimensions and testing date from .csv file    
fn2=fullfile(folder,strcat(extractBefore(basename{1},"_Runtime"),".csv"));
if ~isfile(fn2)
    fprintf('Error: \nCannot find this csv file: %s\n',...
        strcat(extractBefore(basename{1},"_Runtime"),".csv!"))
    return
end
if ~isempty(dir(strcat(folder,extractBefore(basename{1},"_Runtime"),"_*Parsed.xlsx")))
    fprintf('This file %s has been parsed already!\n',extractBefore(basename{1},"_Runtime"))
    return
end
Dimfile=readcell(fn2);
if ismissing(Dimfile{11,8})%file format are different for FAT5,FAT6
    Pos=find(Dimfile{11,3}==',',2,'last');
    widthstr=extractBetween(Dimfile{11,3},Pos(1)+1,Pos(2)-1);
    width=str2double([widthstr{1}]);
    Pos2=find(Dimfile{12,3}==',',2,'last');
    heightstr=extractBetween(Dimfile{12,3},Pos2(1)+1,Pos2(2)-1);
    height=str2double([heightstr{1}]);
else
    width=Dimfile{11,8};
    height=Dimfile{12,8};
end
monthNames = {'January','February','March','April','May','June','July','August','September','October','November','December'}';
%rowDate=find(strcmp(Dimfile(:,2),'date'));
monthnum=find(contains(monthNames,Dimfile{30,5}));
testTime=datetime(str2double(Dimfile{30,7}),monthnum,str2double(Dimfile{30,6}),'F','MM/dd/uuuu');
%TargetStrain=str2double(extractAfter(Dimfile{rowTemp-1,3},","))*2/1e6;
for m=1:size(fullFileName,2)%different temperatures and frequencies for one mixture
    %% fit to sinusoidal function a+c*sin(wt+d)
    Cycle={[]};
    MaximumLoad={[]};
    MinimumLoad={{}};
    Maxload={[]};%direct measurement of peak-to-peak
    Minload={[]};
    MaximumLVDT={[]};
    MinimumLVDT={[]};
    Maxlvdt={[]};%direct peak-to-peak
    Minlvdt={[]};
    p2pstress={[]};
    p2pstrain={[]};
    direct_stress={[]};
    direct_strain={[]};
    flexuralstiff={[]};
    directstiff={[]};%from direct peak-to-peak
    Phaseangle={[]};%from fitting
    Phaseangle2={[]};%hysteresis loop
    %Phaseangle3={[]};% hysteresis loop of fitted curve
    directphaseangle={[]};
    TempC={[]};
    avg_phaseangle={[]};
    avg_E={[]};
    fn=fullFileName{m};
    Specimen=regexp(basename{m},".*\(FAT \d)","match");
    %TargetTemp=regexp(Specimen, '(\d+)C', 'tokens');
    %TargetTemp=str2num(TargetTemp{:});
    TargetFreq=str2double(regexp(basename{m},'(?<=\()[^)]*(?=\)*Hz)', 'match'));

    Rawdata=readcell(fn);
    
    Variablenames=["Cycle Number","Time (sec)","Acurator (mm)","Load (kN)",...
        "Rear  LVDT (mm)","Core Temperature (°C)"];

    %% Look for the variable columns
    
     if find(ismember(Rawdata(1,:),'On-Specimen (mm)'))%fat5 has the rear LVDT as "On-Specimen (mm)" 
         idLVDT=find(ismember(Rawdata(1,:),'On-Specimen (mm)'));
     else 
         idLVDT=find(ismember(Rawdata(1,:),'Rear LVDT (mm)'));
    
     end
     if find(ismember(Rawdata(1,:),'Core Temperature (°C)'))
         idTempC=find(ismember(Rawdata(1,:),'Core Temperature (°C)'));
     end 

    Rawdata2=cell2table(Rawdata(2:end,([1:4,idLVDT,idTempC])),"VariableNames",Variablenames); 
    Rawdata2.("Cycle Number")=fillmissing(Rawdata2.("Cycle Number"),'previous');
    Rawdata2.("Core Temperature (°C)")=fillmissing(Rawdata2.("Core Temperature (°C)"),'previous');    

    n=length(unique(Rawdata2.("Cycle Number")));
    cycnumber=unique(Rawdata2.("Cycle Number"));
    frequency{m}=(Rawdata2.("Cycle Number")(end)-Rawdata2.("Cycle Number")(1))/...
        (Rawdata2.("Time (sec)")(end)-Rawdata2.("Time (sec)")(1));    %obtain the frequency and temperature information
    Temperature{m}=Rawdata2.("Core Temperature (°C)")(1);
    %% Fitting each cycle    
    for i=1:n%n
        cycle=cycnumber(i);
        Time=Rawdata2.("Time (sec)")(Rawdata2.("Cycle Number")==cycle);
        Displacement=Rawdata2.("Rear  LVDT (mm)")(Rawdata2.("Cycle Number")==cycle);
        Load=Rawdata2.("Load (kN)")(Rawdata2.("Cycle Number")==cycle);     
        
        [SineD]=sineFitSSE(Time,Displacement,TargetFreq);%sineFitSSE(fittedx,fittedy,frequency,plots?), 0.05 Hz or 10 Hz
        DisplacementPeakValues=[(pi/2-SineD(3))/(2*pi*SineD(4)),SineD(1)+SineD(2),...
            (pi*3/2-SineD(3))/(2*pi*SineD(4)),SineD(1)-SineD(2)];%[DisplacementPeakTime,DisplacementPeak,DisplacementPeakTime2,DisplacementPeak2]
        %[DisplacementPeakValues]=FatPhaseAngle(Time,SineD);%[DisplacementPeakTime,DisplacementPeak,DisplacementPeakTime2,DisplacementPeak2]
        [SineL]=sineFitSSE(Time,Load,TargetFreq);%sineFitSSE,0.05 Hz or 10Hz
        LoadPeakValues=[(pi/2-SineL(3))/(2*pi*SineL(4)),SineL(1)+SineL(2),...
            (pi*3/2-SineL(3))/(2*pi*SineL(4)),SineL(1)-SineL(2)];%[LoadmaxTime,Loadmax,LoadminTime,Loadmin]
        %[LoadPeakValues] = FatPhaseAngle(Time,SineL);%[LoadmaxTime,Loadmax,LoadminTime,Loadmin]
        %Phaseangle{i}=360*SineD(4)*abs(LoadPeakValues(1)-DisplacementPeakValues(1));
        if SineD(2)>0
            nD=ceil((Time(1)*2*pi*(SineD(4))+SineD(3)-pi/2)/(2*pi));
            TimeD=(pi/2+nD*2*pi-SineD(3))/(2*pi*SineD(4));
        else
            nD=ceil((Time(1)*2*pi*(SineD(4))+SineD(3)+pi/2)/(2*pi));
            TimeD=(-pi/2+nD*2*pi-SineD(3))/(2*pi*SineD(4));
        end
        if SineL(2)>0
            nL=ceil((Time(1)*2*pi*(SineL(4))+SineL(3)-pi/2)/(2*pi));
            TimeL=(pi/2+nL*2*pi-SineL(3))/(2*pi*SineL(4));
        else
            nL=ceil((Time(1)*2*pi*(SineL(4))+SineL(3)+pi/2)/(2*pi));
            TimeL=(-pi/2+nL*2*pi-SineL(3))/(2*pi*SineL(4));
        end
       
        %TempPhase=round((SineL(3)-SineD(3))/(2*pi),0);
        Phaseangle{i}=360*SineD(4)*abs(TimeD-TimeL);
        Displacement2=[Displacement' Displacement(1)]';%connect first and last dot
        
        Load2=[Load' Load(1)]';
        Area2=trapz(Displacement2,Load2);
        Dis0=1/2*(max(Displacement)-min(Displacement));
        Load0=1/2*(max(Load)-min(Load));   
        Phaseangle2{i}=asind(Area2/(pi*Dis0*Load0));
        
        Cycle{i}=cycle;
        MaximumLoad{i}=LoadPeakValues(2);
        [Maxload{i},idx11]=max(Load);
        MinimumLoad{i}=LoadPeakValues(4);
        [Minload{i},~]=min(Load);
        MaximumLVDT{i}=DisplacementPeakValues(2);
        [Maxlvdt{i},idx22]=max(Displacement);
        MinimumLVDT{i}=DisplacementPeakValues(4);
        [Minlvdt{i},~]=min(Displacement);
        Peaktopeakload=abs(MaximumLoad{i}-MinimumLoad{i});
        Peaktopeaklvdt=abs(MaximumLVDT{i}-MinimumLVDT{i});
        direct_peaktopeakload=abs(Maxload{i}-Minload{i});
        direct_peaktopeaklvdt=abs(Maxlvdt{i}-Minlvdt{i});
        p2pstress{i}=0.3555*Peaktopeakload*1000/(width/1000*(height/1000)^2);%unit:Pa
        p2pstrain{i}=12*Peaktopeaklvdt/1000*height/1000/(3*0.3555^2-4*0.1185^2);%unit:m/m (0.3555,0.1185)
        direct_stress{i}=0.3555*direct_peaktopeakload*1000/(width/1000*(height/1000)^2);%unit:Pa
        direct_strain{i}=12*direct_peaktopeaklvdt/1000*height/1000/(3*0.3555^2-4*0.1185^2);%unit:m/m (0.3555,0.1185)
        flexuralstiff{i}=p2pstress{i}*10^-6/p2pstrain{i};%unit:MPa
        directstiff{i}=direct_stress{i}*10^-6/direct_strain{i};
        %TT=Rawdata2.("Core Temperature (°C)")(Rawdata2.("Cycle Number")==cycle);
        %TempC{i}=TT(1);
        avg_E{i}=(flexuralstiff{i}+directstiff{i})/2;
        avg_phaseangle{i}=(Phaseangle{i}+Phaseangle2{i})/2;
        lagtime2=abs(Time(idx11)-Time(idx22));
        directphaseangle{i}=360*10*lagtime2;
        avg_p2pstress{i}=(p2pstress{i}+direct_stress{i})/2e6;%Pa to MPa
        avg_p2pstrain{i}=(p2pstrain{i}+direct_strain{i})/2;
    end
    p2p_stress{m}=mean(cell2mat(avg_p2pstress));%MPa
    p2p_strain{m}=mean(cell2mat(avg_p2pstrain));%mm/mm
    EG{m}=mean(cell2mat(avg_E));%MPa
    PhaseAngle{m}=mean(cell2mat(avg_phaseangle));

end

Frequency=frequency';
P2PStress=p2p_stress';%MPa
P2PStrain=p2p_strain';
E_MPa=EG';%fit sinusoidal model 
Phase_Angle=PhaseAngle';
Temp_C=Temperature';

figure
scatter([Frequency{:}],[E_MPa{:}])

OutputDataSummary={Specimen{1},strcat("FAT",num2str(MachineID)),testTime}';
OLTSfilename="4PB FS ProcessorforOLTS.xlsx";
ProcessedName=strcat(folder,Specimen,"FS_Parsed.xlsx");
copyfile(OLTSfilename,ProcessedName)
JJ=table(Frequency,P2PStress,P2PStrain,E_MPa,Phase_Angle,Temp_C);
writetable(JJ,ProcessedName,'Sheet','Data.Clean','Range','B2','WriteVariableNames',false)
writecell(OutputDataSummary,ProcessedName,'Sheet','Data.Summary','Range','B2');
close all
