function [] = ParseFatUCD(hz)
%To run: ParseFatUCD
% Unless a frequency different from 10Hz:ParseFatUCD(0.05) FOR 0.05 Hz
% require the two output excel files from the UTS015
% require the empty excel processor "4PB Fatigue ProcessorforOLTS.xlsx"
% require function sinFitSSE

if nargin<1
    freq=10;%default frequency value 10 Hz
else 
    freq=hz;
end
[basename,folder]=uigetfile(('*/*RuntimeData*.csv'),'select all the input data files (s)','Multiselect','on');%read raw data file
fullFileName=fullfile(folder,basename);
if ~iscell(fullFileName)
    fullFileName={fullFileName};
    basename={basename};
end
hw = waitbar(0,'Running...');

for m=1:size(fullFileName,2)
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
    Specimen=extractBefore(basename{m},"_Runtime");
    if ~isempty(dir(strcat(folder,Specimen,"_*Parsed.xlsx")))
        fprintf('This file %s has been parsed already!\n',Specimen)
        continue
    end
    Rawdata=readcell(fn);
    MachineID=5;
    if contains(lower(basename{m}),"fat6")
        MachineID=6;
    end
    Variablenames=["Cycle Number","Time (sec)","Acurator (mm)","Load (kN)",...
        "Rear  LVDT (mm)","Core Temperature (°C)"];

    %% Look for the variable columns
    if MachineID==5 % 
         if find(ismember(Rawdata(1,:),'Front LVDT (?)'))%fat5 has the rear LVDT as "front LVDT" or "Rear LVDT" 
             idLVDT=find(ismember(Rawdata(1,:),'Front LVDT (?)'));
         else 
             idLVDT=find(ismember(Rawdata(1,:),'Rear LVDT (mm)'));
    
         end
         if find(ismember(Rawdata(1,:),'Core Temperature (°C)'))
             idTempC=find(ismember(Rawdata(1,:),'Core Temperature (°C)'));
         else 
             idTempC=find(ismember(Rawdata(1,:),'Temperature (°C)'));
    
         end
        else%Fat6
         idLVDT = find(ismember(Rawdata(1,:),'Centre Top LVDT (mm)')|...
             ismember(Rawdata(1,:),'Rear  LVDT (mm)'));%fat6 has the rear LVDT as "Centre Top LVDT" 
         idLVDT=idLVDT(1);
         idTempC= find(ismember(Rawdata(1,:),'Core Temperature (°C)'));
    end
    Rawdata2=cell2table(Rawdata(2:end,([1:4,idLVDT,idTempC])),"VariableNames",Variablenames);
 
  
    Rawdata2.("Cycle Number")=fillmissing(Rawdata2.("Cycle Number"),'previous');
    Rawdata2.("Core Temperature (°C)")=fillmissing(Rawdata2.("Core Temperature (°C)"),'previous');

    % specimen dimensions and testing date from .csv file    
    fn2=fullfile(folder,strcat(Specimen,".csv"));
    if ~isfile(fn2)
        fprintf('Error: \nCannot find this csv file: %s\n',strcat(Specimen,".csv!"))
        continue
    end
    Dimfile=readcell(fn2);
    Pos=find(Dimfile{11,2}==',',2,'last');
    widthstr=extractBetween(Dimfile{11,2},Pos(1)+1,Pos(2)-1);
    width=str2double([widthstr{1}]);
    Pos2=find(Dimfile{12,2}==',',2,'last');
    heightstr=extractBetween(Dimfile{12,2},Pos2(1)+1,Pos2(2)-1);
    height=str2double([heightstr{1}]);
    monthNames = {'January','February','March','April','May','June','July','August','September','October','November','December'}';
    rowDate=find(strcmp(Dimfile(:,2),'date'));
    monthnum=find(contains(monthNames,Dimfile{rowDate,5}));
    testTime=datetime(str2double(Dimfile{rowDate,7}),monthnum,str2double(Dimfile{rowDate,6}),'F','MM/dd/uuuu');
    rowTemp=find(strcmp(Dimfile(:,2),'temperature'));
    TargetTemp=str2double(extractAfter(Dimfile{rowTemp,3},","));
    TargetFreq=str2double(extractAfter(Dimfile{rowTemp+1,2},","));
    TargetStrain=str2double(extractAfter(Dimfile{rowTemp+3,6},","))/1e6;


    %% obtain delayed time from first 10 cycles fitting
    DelayData=Rawdata2(Rawdata2.("Cycle Number")<11,:);
    [DelayFit]=sineFitSSE(DelayData.("Time (sec)"),DelayData.("Load (kN)"),freq);
    DelaySec=(1/DelayFit(4)-DelayFit(3)/(2*pi*DelayFit(4)));
    DelayPoint=size(DelayData(DelayData.("Time (sec)")<DelaySec,1),1);
    totrow=size(Rawdata2(Rawdata2.("Time (sec)")>=DelaySec,1),1);
    totcyc=Rawdata2.("Cycle Number")(1:totrow);
    Rawdata2=Rawdata2(Rawdata2.("Time (sec)")>=DelaySec,:);
    Rawdata2.("Cycle Number")=totcyc;
    
    Rawdata2=Rawdata2(Rawdata2.("Cycle Number")<(totcyc(end)-5),:);%delete last five cycles due to unstable sino reading
    n=length(unique(Rawdata2.("Cycle Number")));
    cycnumber=unique(Rawdata2.("Cycle Number"));
    %% Fitting each cycle    
    for i=1:n%n
        cycle=cycnumber(i);
        Time=Rawdata2.("Time (sec)")(Rawdata2.("Cycle Number")==cycle);
        Displacement=Rawdata2.("Rear  LVDT (mm)")(Rawdata2.("Cycle Number")==cycle);
        Load=Rawdata2.("Load (kN)")(Rawdata2.("Cycle Number")==cycle);

        
%         if i>1 && i<n && (cycnumber(i+1)-cycle)>1 %sequence jump
%             Time=Time(1:(end-DelayPoint));
%             Displacement=Displacement(1:(end-DelayPoint));
%             Load=Load(1:(end-DelayPoint));
%         end         

        [SineD]=sineFitSSE(Time,Displacement,freq,0);%sineFitSSE(fittedx,fittedy,frequency,plots?), 0.05 Hz or 10 Hz
        DisplacementPeakValues=[(pi/2-SineD(3))/(2*pi*SineD(4)),SineD(1)+SineD(2),...
            (pi*3/2-SineD(3))/(2*pi*SineD(4)),SineD(1)-SineD(2)];%[DisplacementPeakTime,DisplacementPeak,DisplacementPeakTime2,DisplacementPeak2]
        %[DisplacementPeakValues]=FatPhaseAngle(Time,SineD);%[DisplacementPeakTime,DisplacementPeak,DisplacementPeakTime2,DisplacementPeak2]
        [SineL]=sineFitSSE(Time,Load,freq,0);%sineFitSSE,0.05 Hz or 10Hz
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
        TT=Rawdata2.("Core Temperature (°C)")(Rawdata2.("Cycle Number")==cycle);
        TempC{i}=TT(1);
        avg_E{i}=(flexuralstiff{i}+directstiff{i})/2;
        avg_phaseangle{i}=(Phaseangle{i}+Phaseangle2{i})/2;
        lagtime2=abs(Time(idx11)-Time(idx22));
        directphaseangle{i}=360*10*lagtime2;
        p2pstress{i}=p2pstress{i}/1e6;%Pa to MPa

    
    end
    Cycle_Number=Cycle';
%     MaxLoad=MaximumLoad';
%     DirectMaxLoad=Maxload';
%     MinLoad=MinimumLoad';
%     DirectMinLoad=Minload';
%     MaxLVDT=MaximumLVDT';
%     DirectMaxLVDT=Maxlvdt';
%     MinLVDT=MinimumLVDT';
%     DirectMinLVDT=Minlvdt';
    p2p_stress=p2pstress';%MPa
%    Directstress=direct_stress';
    p2p_strain=p2pstrain';
%    Directstrain=direct_strain';
    E_MPa=flexuralstiff';%fit sinusoidal model 
%     Direct_E=directstiff';% UTS 015 output
%     Phase_Angle=Phaseangle';
%     Phase_Angle2=Phaseangle2';
%     Direct_phaseangle=directphaseangle';
    Temp_C=TempC';
    avg_angle=avg_phaseangle';
    avg_stiff=avg_E';
    figure
    scatter([Cycle_Number{:}],[E_MPa{:}])
    hold on
    scatter([Cycle_Number{:}],[avg_stiff{:}])
    if mod(i,floor(n/10))<1e-2
     waitbar((m/size(fullFileName,2)),hw);
    %waitbar(((m-1)*n+i)/(n*size(fullFileName,2)),hw);
    end

    
%     JJ=table(Cycle_Number,p2p_stress,p2p_strain,E_MPa,Phase_Angle,Phase_Angle2,Directstress,...
%        Directstrain,Direct_E,Direct_phaseangle,MaxLoad,MinLoad,...
%        MaxLVDT,MinLVDT,DirectMaxLoad,DirectMinLoad,DirectMaxLVDT,...
%        DirectMinLVDT,Temp_C,avg_stiff,avg_angle);
%     JJ.Properties.VariableNames=["cycle","p2pStress(Pa)","p2pstrain (m/m)"...
%         "E/G(MPa)","Phase Angle","Phase Angle2","DirectStress(Pa)","DirectStrain(m/m)","DirectE(MPa)",...
%         "Direct_PhaseAngle","MaxLoad","MinLoad","MaxLVDT","MinLVDT","DirectMaxLoad",...
%         "DirectMinLoad","DirectMaxLVDT","DirectMinLVDT","Temp(C)","Avg.E","Avg.angle"];
%     writetable(JJ,strcat(Specimen,"_processed2.xlsx"),"FileType",...
%         "spreadsheet", 'WriteVariableNames',true)

%% calculate variables: E50, Date,Nf_50,Nf_20,Nf_nxSR for Data.Summary
  TestMachine=strcat("FAT",num2str(MachineID));
  E50=avg_E{[Cycle{:}]==50};
  Temp=find(([avg_E{:}]/E50-0.5)<0);
  Nf50=Cycle{Temp(1)};
  Temp2=find(([avg_E{:}]/E50-0.2)<0);%extrapolation for 20%SR
  if isempty(Temp2) %some tests used the 15% cycleXstiff
    SR20cycle=[Cycle{(end-11):(end-1)}];
    SR20SR=[avg_E{(end-11):(end-1)}]/E50;  
  else
    SR20cycle=[Cycle{(Temp2(1)-11):(Temp2(1)-1)}];
    SR20SR=[avg_E{(Temp2(1)-11):(Temp2(1)-1)}]/E50;  
  end
  Nf20=interp1(SR20SR,SR20cycle,0.2,'linear','extrap');
  nxSR=[Cycle{:}].*[avg_E{:}];
  NfnxSR=Cycle{nxSR==max(nxSR)};
  OutputDataSummary={Specimen,TestMachine,testTime,E50,Nf50,NfnxSR,Nf20,...
     TargetTemp,TargetFreq,TargetStrain,TempC{1},p2pstrain{50}}';
%% write into 4PB processor for OLTS
  OLTSfilename="4PB Fatigue ProcessorforOLTS.xlsx";
  ProcessedName=strcat(folder,Specimen,"_",num2str(round(TargetStrain*1e6)),"ue_Parsed.xlsx");
  copyfile(OLTSfilename,ProcessedName)
  writecell(Cycle_Number,ProcessedName,'Sheet','Data.Clean','Range','A2');
  JJ=table(p2p_stress,p2p_strain,avg_stiff,avg_angle,Temp_C);
  writetable(JJ,ProcessedName,'Sheet','Data.Clean','Range','C2','WriteVariableNames',false)
  writecell(OutputDataSummary,ProcessedName,'Sheet','Data.Summary','Range','B2');
  close all
 
end
  delete(hw);
  %load handel
  %sound(y,Fs)
end