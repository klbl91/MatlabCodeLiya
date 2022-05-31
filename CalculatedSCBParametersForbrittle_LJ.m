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
%   Gauss1 is used here for fitting and the inflection point is the local minimum 
%   of second derivative of Gauss1 fitting result 
%------------------------------------------------------------------

%% calculation steps need to update for Index=1 and Index=0 following
%%Index=2 part
clc;
clear;
%cd 'D:\1.study\PHD\UCPRC\Project\SCB\SCB Test data\SCB_test data archive\ILL\ANALYSIS CODES AND RESULTS';%% to update
%%
%%input data file, information about the geometry
%testdatapath="D:\1.study\PHD\UCPRC\Project\SCB\SCB Test data\464_LMLC_test_data";%% to update
testdatapath=pwd;
Index=0;%Index ==0 calculate the parameters requiring geometry information, otherwise Index=1 OR 2 depending on which parameters to obtain
%%
if Index==0
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
           error(msg);
        end
    end
    disp("The specimens that need to be calculated are:")
    disp(s.textdata(:,k));

    %%read the SCB test data
%     disp('Please press a key to continue')  
%     pause;

    [basename2,folder2]=uigetfile(testdatapath+"\*.csv",'select the SCB testing results file (s)','Multiselect','on');
    fullFileName2=fullfile(folder2,basename2);
    TestsampleID={[]};
    Gf={[]};
    FI={[]};
    FIasc={[]};
    Spp={[]};
    Sasc={[]};
    frontarea={[]};
    backarea={[]};
    Strength={[]};
    KIC={[]};
    if size(fullFileName2,2)~=size(s.data,1)
        msg="The testing files number does not match with the input file samples number";
        error(msg);
    end

    for m=1:size(fullFileName2,2)
        DATA=xlsread(fullFileName2{m});
        TestsampleID{m}=extractBefore(basename2{m},' ');
        [peakload,peakrow]=max(DATA(:,3));
        Strength{m}=peakload*10^9/(2*s.data(m,r)*s.data(m,j));%%unit N/m2
        %%fit the testing result curves;
        %%before the peak: a polynomial equation with a degree of six
        %%P1(u)=C1*u^6+C2*u^5+C3*u^4+C4*u^3+C5*u^2+C6*u+C7;
        %%after the peak: a exponential equation with 4 terms
        %%P2(u)=sum(di*exp[-((u-ei)/fi)^2]
        preU=DATA(1:peakrow,4);
        preLoad=DATA(1:peakrow,3);
        postU=DATA(peakrow:length(DATA),4);
        postLoad=DATA(peakrow:length(DATA),3);
        frontarea{m}=trapz(preU,preLoad);
        backarea{m}=trapz(postU,postLoad);
        
       
        %%calculate Gf
        p1=polyfit(preU,preLoad,6);
        q1=polyint(p1);
        Wf1=diff(polyval(q1,[0,DATA(peakrow,4)]));

        p2=fit(postU,postLoad,'gauss2');
        %DispEnd=arrayfun(@(y)fzero(@(x)p2(x)-y,1),0.1);%calculate the corresponding displacement when loading =0.1kN
        [~,INDEXdispEnd]=min(abs(DATA(peakrow:length(DATA),3)-0.1));
        DispEnd=DATA(INDEXdispEnd+peakrow-1,4);
        Wf2=abs(integrate(p2,DATA(peakrow,4),DispEnd));
        Wf=Wf1+Wf2;
        Gf{m}=Wf*10^6/((s.data(m,i))*s.data(m,j));
        %%calculate slope
        DATA1=DATA(1:peakrow,:);
        [~,point1r]=min(abs(DATA1(:,3)-1/4*DATA1(peakrow,3)));
        [~,point2r]=min(abs(DATA1(:,3)-3/4*DATA1(peakrow,3)));
        Sasc{m}=(DATA1(point2r,3)-DATA1(point1r,3))/(DATA1(point2r,4)-DATA1(point1r,4));
        syms x
        eq2=p2.a1*exp(-((x-p2.b1)/p2.c1)^2);
        eq3=diff(eq2);
        p22=diff(diff(eq2));%second derivative
        PeakDisp=DATA(peakrow,4);
        tt=1;
        for incr=PeakDisp:0.0001:DispEnd
            inflect(tt)=abs(double(subs(p22,x,incr)));
            tt=tt+1;
        end
        localMin=islocalmin(inflect);
        inflectIndex=find(localMin==1);
        inflect_pt=PeakDisp+0.0001*inflectIndex;
        [~,inflect_ptindex]=min(abs(postU-inflect_pt));
        inflect_load=postLoad(inflect_ptindex);
        Spp{m}=double(subs(eq3,x,inflect_pt));
        %%calculate the FI
        FI{m}=Gf{m}*0.01/abs(Spp{m});
        FIasc{m}=Gf{m}*0.01/Sasc{m};
        %%calculate KIC
        YI=4.782+1.219*(s.data(m,r)-s.data(m,i))/s.data(m,r)+0.063*exp(7.045*(s.data(m,r)-s.data(m,i))/s.data(m,r));    %TT.YI=4.782+1.219*(TT.Notch./TT.Radius)+0.063*exp(7.045*(TT.Notch./TT.Radius));
    %TT.KICmatlab=TT.YI.*sqrt(pi()*TT.Notch/1000).*(TT.Pmax*1e3./(2*TT.Radius.*TT.TT.Thickness);Notch,Radius and thickness are in mm; Pmax in kN                                                                    

        KIC{m}=YI*sqrt(pi()*(s.data(m,r)-s.data(m,i))/1000)*peakload*10^3/(2*s.data(m,r)*s.data(m,j));%the unit of KIC is MPa(sqrt(m))
                
           
        figure
        plot(DATA(:,4),DATA(:,3));
        title(TestsampleID{m});
        xlabel('Deformation (mm)');
        ylabel('Load (kN)');
        hold on
        syms xx
        y = Spp{m}*(xx-inflect_pt)+inflect_load;
        fplot(y, [PeakDisp DispEnd], '--r')
        plot(inflect_pt,inflect_load,'r*')
        ylim([0 peakload])
        hold off
        saveas(gcf,strcat(TestsampleID{m},'.png'))
        
    end
     %%Fracture engery, FI, FIasc, front slope,
    %%       back slope, KIC, strength, prepeak area, afterpeak area)

    SampleID=TestsampleID';
    FractureEnergy=Gf';
    FImatlab=FI';
    FIascmatlab=FIasc';
    Sppmatlab=Spp';
    Sascmatlab=Sasc';
    KICmatlab=KIC';
    strength=Strength';
    Area1=frontarea';
    Area2=backarea';
    JJ=table(SampleID,FractureEnergy,FImatlab,FIascmatlab,Sppmatlab,Sascmatlab,...
        KICmatlab,...
        strength,Area1,Area2);
    JJ.Properties.VariableNames={'Sample_ID','Fracture_energy','FImatlab','FIascmatlab',...
        'Sppmatlab','Sascmatlab','KICmatlab','strength','Area1','Area2'};
    JJ.Properties.VariableUnits={' ','J/m^2',' ',' ',' ',' ','MPa/m^1/2','N/m^2',...
        'N*m','N*m'};
    CellID2=char('A'+size(s.data,2));
    StartCell2=strcat(CellID2,num2str(1));
    writetable(JJ,fullFileName,'FileType','spreadsheet','Range',StartCell2);
end
if Index==1
    [basename,folder]=uigetfile(testdatapath+"\*.csv",'select the test data files','Multiselect','on');
    fullFileName=fullfile(folder,basename);
    for m=1:size(fullFileName,2)
        DATA=xlsread(fullFileName{m});
        TestsampleID{m}=extractBefore(basename,' ');
        LoadTimesDisp=DATA(:,3).*DATA(:,4);
        [MAX,INDEX]=max(LoadTimesDisp);
        MaxEnergy{m}=MAX;
        MaxDisp{m}=DATA(INDEX,4);%the displacement when the peak of load x displacement accurs.
    end
       SpecimenID=TestsampleID{1,1}';
       LoadTimesDispPeak=MaxEnergy';
       DispAtPeak=MaxDisp';
       JJ=table(SpecimenID, LoadTimesDispPeak,DispAtPeak);
       writetable(JJ,'New parameters.xlsx','FileType','spreadsheet','Sheet','GLE5-FL');
    
end

if Index==2
    [basename,folder]=uigetfile(testdatapath+"\*.csv",'select the test data files','Multiselect','on');
    h=waitbar(0,'Please wait...');
    fullFileName=fullfile(folder,basename);
    for m=1:size(fullFileName,2)
        DATA=xlsread(fullFileName{m});
        TestsampleID{m}=extractBefore(basename,' ');
        [PeakLoad,index]=max(DATA(:,3));
        PeakDisp=DATA(index,4);
        PostPeakdata=DATA(index:length(DATA),:);
        [~,index75]=min(abs(PostPeakdata(:,3)-0.75*PeakLoad));
        L75{m}=PostPeakdata(index75,4);%the displacement of L75 (corresponding displacement at 75% of peak load at the post curve)
        %calculate Spp and related displacement
        postU=DATA(index:length(DATA),4);
        postLoad=DATA(index:length(DATA),3);
        [~,INDEXdispEnd]=min(abs(DATA(index:length(DATA),3)-0.1));
        DispEnd=DATA(INDEXdispEnd+index-1,4);
%         flowbound=[0.9,0.7,0.5,0.3,0.1,0.05,0.01,0.005,0.001];
%         for bound=1:9
%         f1=flowbound(bound);   
%         options=fitoptions('gauss4','Lower',[-10 -10 f1 -10 -10 f1 -10 -10 f1...
%         -10 -10 f1],'Upper',[10 10 10 10 10 10 10 10 10 10 10 10]);   
%         [p2,pof]=fit(postU,postLoad,'gauss4',options);
%         Squareroot(bound)=pof.sse;
%         end
%         [~,index]=min(Squareroot);
%         f2=flowbound(index);   
%         options=fitoptions('gauss4','Lower',[-10 -10 f2 -10 -10 f2 -10 -10 f2...
%         -10 -10 f2],'Upper',[10 10 10 10 10 10 10 10 10 10 10 10]);   
%         p2=fit(postU,postLoad,'gauss4',options);
%         syms x
%         eq2=p2.a1*exp(-((x-p2.b1)/p2.c1)^2)+p2.a2*exp(-((x-p2.b2)/p2.c2)^2)+...
%         p2.a3*exp(-((x-p2.b3)/p2.c3)^2)+p2.a4*exp(-((x-p2.b4)/p2.c4)^2);
        [xData, yData] = prepareCurveData( postU, postLoad );        
        ft = fittype( 'gauss4' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
        %opts.StartPoint = [3.32281 1.07718 0.160188647019638 2.76343756258727 1.36976 0.15395657885066 2.01229955245587 1.63615 0.163969669505336 1.39628447338922 1.91561 0.185652310144618];
        opts.StartPoint = [3 1 0.2 3 1.5 0.15 2.5 2 0.15 1.7 2.3 0.15];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        syms x
        eq2=fitresult.a1*exp(-((x-fitresult.b1)/fitresult.c1)^2)+fitresult.a2*exp(-((x-fitresult.b2)/fitresult.c2)^2)+...
         fitresult.a3*exp(-((x-fitresult.b3)/fitresult.c3)^2)+fitresult.a4*exp(-((x-fitresult.b4)/fitresult.c4)^2);
        eq3=diff(eq2);
        p22=diff(diff(eq2));%second derivative
        gg2=matlabFunction(p22);
        inflect_pt{m}=fzero(gg2,PeakDisp);%displacement at inflection point
        inflect_load{m}=double(subs(eq2,x,inflect_pt{m}));       
        Spp{m}=double(subs(eq3,x,inflect_pt{m}));
        CrossDisplace{m}=inflect_pt{m}+inflect_load{m}/abs(Spp{m});
        plot(DATA(:,4),DATA(:,3))
        hold on
        fplot(eq2,[PeakDisp,DispEnd])
        hold on
        for iii=1:100
            if Spp{m}>0|CrossDisplace{m}>DispEnd|inflect_pt{m}<PeakDisp
                inflect_pt{m}=fzero(gg2,PeakDisp+iii*0.1);%0.1,displacement at inflection point
                inflect_load{m}=double(subs(eq2,x,inflect_pt{m}));       
                Spp{m}=double(subs(eq3,x,inflect_pt{m}));
                CrossDisplace{m}=inflect_pt{m}+inflect_load{m}/abs(Spp{m});
            else
                break
            end
        end
        CrossDisplace{m}=inflect_pt{m}+inflect_load{m}/abs(Spp{m});
        plot(inflect_pt{m},inflect_load{m},'ko')
        hold on
        xxx=PeakDisp:0.2:DispEnd;
        yyy=Spp{m}*(xxx-inflect_pt{m})+inflect_load{m};
        plot(xxx,yyy)
        xlim([0 DispEnd])
        ylim([0 PeakLoad])
        hold off

        waitbar(m*100/size(fullFileName,2),h)
    end
       close(h)
       SpecimenID=TestsampleID{1,1}';
       L=L75';
       inflect_displace=inflect_pt';
       SPP=Spp';
       inflect_force=inflect_load';
       ExtendDisplace=CrossDisplace';
       JJ=table(SpecimenID, SPP,L,inflect_displace,inflect_force,ExtendDisplace);
       JJ.Properties.VariableNames={'Sample_ID','Spp','L75','Displacement_inflectionPoint',...
           'Load_inflectionpoint','Displacement_extendtoXaxis'};
       writetable(JJ,'SCB_Displacements.xlsx','FileType','spreadsheet','Sheet','464SCBdisplacements3');
    
end
