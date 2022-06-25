%% This module is to used for data reduction
clc;clearvars -except Final_Candidates
Heptacosa_fragments=[68.99521; 168.988822; 230.991896; 263.987106; 280.988703; 301.983912; 382.979122; 401.977525; 413.977525; 451.974332; 475.974332; 482.9733; 513.971138; 594.9669; 613.9653; 632.963154];
% load('Final_Candidates.mat');
T0=6;                         % scan number of RT=6 min
T_end=51.5;                   % scan number of RT=51.5 min
load('RUSBoosted_APGC_QToF.mat');
load('ID_library.mat')
load('XIC_Primary.mat')
ID_I=cell2mat(XIC_Primary(:,4));
ID_E=cell2mat(XIC_Primary(:,5));
%%
Delta_mass=0.01;  % (Da)
% Delta_RT=0.1; % (min)
AAA=zeros(4e5,7);
counter=0;c=0;
% load('C:\Users\fakours\Desktop\APGC\Noise_Neg_mz_RT.mat') % detected features in blank
for i=1:size(Final_Candidates,1)
    c=c+1;
    if c==1000
        i
        c=0;
    end
    if ~ischar(Final_Candidates{i,1})
        ID_M=Final_Candidates{i,1};
        x_N_Iso=find(ID_M>=ID_I & ID_M<=ID_E);
        N_Iso=XIC_Primary{x_N_Iso,2};       % Number of isotopologues
    end
    if isempty(Final_Candidates{i,3}) || Final_Candidates{i,3}==1
        x=find(abs(Heptacosa_fragments-Final_Candidates{i,4})<=Delta_mass);
        % y=find(abs(Noise_Neg_mz_RT(:,1)-Final_Candidates{i,4})<=Delta_mass & abs(Noise_Neg_mz_RT(:,2)-Final_Candidates{i,9})<=Delta_RT);
        if isempty(x) && ... % isempty(y) && ...
                Final_Candidates{i,9}<51 && ...
                Final_Candidates{i,11}>2 && ...
                Final_Candidates{i,13}<0.8 && ...
                Final_Candidates{i,12}>25
            %% Si(CH3)O Mass defect
            M=Final_Candidates{i,4};
            MD=59/58.99533*M-floor(M);
            MD_Cond=1;
            if MD>0.08 && MD<0.2 && M<250
                MD_Cond=0;
            elseif MD+1>1.08 && MD+1<1.2 && M<250
                MD_Cond=0;
            end
            if  MD_Cond==1
                %% Boiling point prediction
                cl=ID_library(ID_M,4);br=ID_library(ID_M,3);
                if cl>0 && br==0
                    M=Final_Candidates{i,4}+35-16;
                    PTb=BP_Predictor(M,cl+1,0);      % Predicted boiling point
                elseif br>0 && cl==0
                    M=Final_Candidates{i,4}+79-16;
                    PTb=BP_Predictor(M,0,br+1);      % Predicted boiling point
                elseif cl>0 && br>0
                    M=Final_Candidates{i,4}+79-16;
                    PTb=BP_Predictor(M,cl,br+1);     % Predicted boiling point
                end
                %% RT boundaries prediction
                PRT1=(0.8*PTb-96.1880)/15.5940; % Predicted the front end of RT
                PRT2=(1.2*PTb+64.5008)/9.25;    % Predicted the end of RT
                if PRT2>Final_Candidates{i,9} && PRT1<Final_Candidates{i,9}
                    %% False positive prediction module
                    S = RUSBoosted_APGC_QToF.predictFcn([N_Iso;       % Number of isotopologues
                        Final_Candidates{i,7};  % NEME
                        Final_Candidates{i,8};  % PCS
                        Final_Candidates{i,11}; % NDCS
                        Final_Candidates{i,12}; % RCS
                        Final_Candidates{i,13}]');    % RPW
                    if S==1
                        counter=counter+1;
                        % AAA=[Number of row in the spreadsheet, Compound ID, Mass, Intensity, Retention time, Peak area, Identification score]
                        AAA(counter,:)=[i+1,ID_M,Final_Candidates{i,4},Final_Candidates{i,6},Final_Candidates{i,9},Final_Candidates{i,10},Final_Candidates{i,14}];
                    end
                end
            end
        end
    end
end
AAA=AAA(1:counter,:);
AAA=sortrows(AAA,1);
%% To remove [M+1], [M+2], [M+3], ... peaks
MplusOne=[0;1.9974;3.0007;3.9948];
Delta_mass=0.01;            % Mass tolerance (Da)
Delta_RT=0.1;               % Retention time tolerance (min)
rt_HGrid=T0:Delta_RT:T_end;  % to create a meshgrid for the RT
E=zeros(4e3,7);
CounterE=0;
for t=2:length(rt_HGrid)
    x_rt=find(AAA(:,5)>=rt_HGrid(t-1) & AAA(:,5)<rt_HGrid(t));
    if ~isempty(x_rt)
        B=AAA(x_rt,:);
        while size(B,1)~=0
            x_MW=[];
            for M=1:length(MplusOne)
                x_MW=[x_MW;find(abs(B(1,3)-MplusOne(M)-B(:,3))<=Delta_mass)];
                x_MW=[x_MW;find(abs(B(1,3)+MplusOne(M)-B(:,3))<=Delta_mass)];
            end
            x_MW=unique(x_MW);
            H_score=max(B(x_MW,7));
            for i=1:length(x_MW)
                if B(x_MW(i),7)>=H_score*1
                    CounterE=CounterE+1;
                    E(CounterE,:)=B(x_MW(i),:);
                end
            end
            B(x_MW,:)=[];
        end
    end
end
E=E(1:CounterE,:);
rt_HGrid=T0+Delta_RT/2:Delta_RT:T_end+Delta_RT/2;
F=zeros(1e3,7);
CounterF=0;
for t=2:length(rt_HGrid)
    x_rt=find(E(:,5)>=rt_HGrid(t-1) & E(:,5)<rt_HGrid(t));
    if ~isempty(x_rt)
        B=E(x_rt,:);
        while size(B,1)~=0
            x_MW=[];
            for M=1:length(MplusOne)
                x_MW=[x_MW;find(abs(B(1,3)-MplusOne(M)-B(:,3))<=Delta_mass)];
                x_MW=[x_MW;find(abs(B(1,3)+MplusOne(M)-B(:,3))<=Delta_mass)];
            end
            x_MW=unique(x_MW);
            H_score=max(B(x_MW,7));
            for i=1:length(x_MW)
                if B(x_MW(i),7)>=H_score*1
                    CounterF=CounterF+1;
                    F(CounterF,:)=B(x_MW(i),:);
                end
            end
            B(x_MW,:)=[];
        end
    end
end
F=F(1:CounterF,:);
%%
F=sortrows(F,-3);
Reduced_data=[];
while size(F,1)~=0
    x_mz_rt=find(abs(F(:,3)-F(1,3))<=Delta_mass & abs(F(:,5)-F(1,5))<=Delta_RT);
    [~,x_S]=max(F(x_mz_rt,7));
    Reduced_data=[Reduced_data;F(x_mz_rt(x_S),:)];
    F(x_mz_rt,:)=[];
end
Reduced_data=sortrows(Reduced_data,-3);
save('Reduced_data.mat','Reduced_data')
%%
Halomap_M_X_O_Neg
