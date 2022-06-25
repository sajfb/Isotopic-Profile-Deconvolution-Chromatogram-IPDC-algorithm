%% This section isolates isomers, calculates peak area, and presents candidate molecular formulas
% This module returns molecular formulas for [M]+ ionization pathways
clc;clearvars -except Final_Matches
format longG
route='C:\Users\fakours\Desktop\APGC\';     % the address of essential files
% load('Final_Matches.mat')
%% Chromatographic analysis
min_INT=500;                % Intensity threshold
Chromatogram_scan_interval=10;   % 10 scan times equals to ~ 2.88 sec
Smoothing_window=10;
Max_No_isomers=25;
Mass_error=15/1e3;         % mass error= 15mDa
coeff=[2.8400234831893485	27.09471679378884	0.46361580527372137	0.4787634170258004	0.08426398342725516	0.5690217702234079];
%%
load('XIC_Combination.mat')
load('ID_library.mat')
load('Peaks_TIC1.mat')
load('Retention_Time_TIC1.mat')
EntireIDs=Final_Matches(:,1);
INT=Final_Matches(:,4);
ScanNumber=Final_Matches(:,7);
NEMD_ave=Final_Matches(:,5);
PCS_ave=Final_Matches(:,6);
load('ID_Mass.mat')
load('XIC_Primary.mat')
ID_I=cell2mat(XIC_Primary(:,4));
ID_E=cell2mat(XIC_Primary(:,5));
%%
Column_xlx=[15 19]; % Columns where the molecular formulas were placed
Final_Candidates=cell(2e6,22);
x_MolF=1;
for ID_i=unique(EntireIDs,'stable')'
    x=find(EntireIDs==ID_i);
    if length(x)>2  % Each candidate compounds should be detected in at least 3 chromatogram scans
        x_N_Iso=find(ID_i>=ID_I & ID_i<=ID_E);
        N_Iso=XIC_Primary{x_N_Iso,2};       % Number of isotopologues
        Q=XIC_Primary{x_N_Iso,1};
        IsotopeModel=[Q(:,1)+ID_Mass(ID_i,2),Q(:,3)];
        if N_Iso>3
            noise_removal=1;
            Max_No_isomers=25;
        else
            noise_removal=10;
            Max_No_isomers=5;
        end
        [No_Isomer, ScanNumber_isomers, PA_isomers, n_w500, P_W_r, neme, pcs]=PeakArea(Peaks_TIC1,Retention_Time_TIC1,IsotopeModel,Mass_error,...
            INT(x),ScanNumber(x),Chromatogram_scan_interval,Smoothing_window,noise_removal,Max_No_isomers,min_INT);
        %%
        Cond=0;
        N_I=0;ik=0;PA=0;ss={};w5=0;P_W=0;neme_1=0;pcs_1=0;
        for i=1:No_Isomer
            if length(ScanNumber_isomers{i})>2
                ik=ik+1;
                N_I=N_I+1;
                PA(ik,1)=PA_isomers(i);
                ss(ik,1)=ScanNumber_isomers(i);
                w5(ik,1)=n_w500(i);
                P_W(ik,1)=P_W_r(i);
                neme_1(ik,1)=neme(i);
                pcs_1(ik,1)=pcs(i);
            end
        end
        No_Isomer=N_I;
        if No_Isomer>0
            PA_isomers=PA;
            ScanNumber_isomers=ss;
            n_w500=w5;
            P_W_r=P_W;
            neme=neme_1;
            pcs=pcs_1;
            if No_Isomer==1
                ScanNumber_isomer=ScanNumber_isomers{1};
                x=zeros(size(ScanNumber_isomer,1),1);
                for k=1:size(ScanNumber_isomer,1)
                    x(k)=find(ScanNumber_isomer(k)==ScanNumber & EntireIDs==ID_i,1);
                end
            end
            Cond=1;
        end
        if Cond==1
            %%
            No_Hist=0;
            Final_Candidates{x_MolF,1}=ID_i;  % ID of the compounds in MATLAB code
            Final_Candidates{x_MolF,2}=['[',Chemical_structure_Print(ID_library(ID_i,:)),']+'];   % Molecular Formula Ion
            Final_Candidates{x_MolF,3}=No_Isomer;
            [~,y]=max(INT(x));
            Final_Candidates{x_MolF,5}=Final_Matches(x(1),3);  % Exact mass
            Final_Candidates{x_MolF,10}=round(sum(PA_isomers),2);     % Total peak area of all isomers
            if No_Isomer==1
                Final_Candidates{x_MolF,4}=Final_Matches(x(y(1)),2);  % Experimental mass
                Final_Candidates{x_MolF,6}=Final_Matches(x(y(1)),4);  % Intensity
                Final_Candidates{x_MolF,7}=round(neme,2);  % NEME (mDa)
                Final_Candidates{x_MolF,8}=round(pcs,2);  % PCS (%)
                Final_Candidates{x_MolF,9}=round(Retention_Time_TIC1(Final_Matches(x(y(1)),7)),3);  % Retention time (Min)
                Final_Candidates{x_MolF,11}=length(cell2mat(ScanNumber_isomers));    % Number of detected chromatogram scan(s) (NDCS)
                Final_Candidates{x_MolF,12}=round(Final_Candidates{x_MolF,11}/n_w500*100,2); % RCS
                Final_Candidates{x_MolF,13}=round(P_W_r,2);  % Ratio of the peak width at healf height to the peak width at baseline
                Final_Candidates{x_MolF,14}=round((N_Iso^coeff(1))*(log10(Final_Candidates{x_MolF,6}/500)^1)*...
                    (log10(Final_Candidates{x_MolF,10}/13.96)^1)*((Final_Candidates{x_MolF,8}/1000)^coeff(2))*...
                    (Final_Candidates{x_MolF,11}^coeff(3))*(Final_Candidates{x_MolF,12}^coeff(4))/...
                    ((Final_Candidates{x_MolF,7}/10)^coeff(5))/(Final_Candidates{x_MolF,13}^coeff(6)),0);
            else    % No_Isomer>1
                % this section is to show each isomer
                for i=1:No_Isomer
                    ScanNumber_isomer=ScanNumber_isomers{i};
                    x_MolF=x_MolF+1;
                    Final_Candidates{x_MolF,1}=[num2str(ID_i),'_',num2str(i)];
                    x=zeros(size(ScanNumber_isomer,1),1);
                    for k=1:size(ScanNumber_isomer,1)
                        x(k)=find(ScanNumber_isomer(k)==ScanNumber & EntireIDs==ID_i,1);
                    end
                    [~,y]=max(INT(x));
                    Final_Candidates{x_MolF,4}=Final_Matches(x(y(1)),2);  % Experimental mass
                    Final_Candidates{x_MolF,5}=Final_Matches(x(y(1)),3);  % Exact mass
                    Final_Candidates{x_MolF,6}=Final_Matches(x(y(1)),4);  % Intensity
                    Final_Candidates{x_MolF,7}=round(neme(i),2);  % NEME (mDa)
                    Final_Candidates{x_MolF,8}=round(pcs(i),2);  % PCS (%)
                    Final_Candidates{x_MolF,9}=round(Retention_Time_TIC1(Final_Matches(x(y(1)),7)),3); % Retention time (Min)
                    Final_Candidates{x_MolF,10}=round(PA_isomers(i),2);   % Indiviual peak area of each isomer
                    Final_Candidates{x_MolF,11}=length(ScanNumber_isomer);    % Number of detected chromatogram scan(s) (NDCS)
                    No_Hist=No_Hist+Final_Candidates{x_MolF,11};
                    Final_Candidates{x_MolF,12}=round(Final_Candidates{x_MolF,11}/n_w500(i)*100,2); % RCS
                    Final_Candidates{x_MolF,13}=round(P_W_r(i),2);  % Ratio of the peak width at healf height to the peak width at baseline
                    Final_Candidates{x_MolF,14}=round((N_Iso^coeff(1))*(log10(Final_Candidates{x_MolF,6}/500)^1)*...
                        (log10(Final_Candidates{x_MolF,10}/13.96)^1)*((Final_Candidates{x_MolF,8}/1000)^coeff(2))*...
                        (Final_Candidates{x_MolF,11}^coeff(3))*(Final_Candidates{x_MolF,12}^coeff(4))/...
                        ((Final_Candidates{x_MolF,7}/10)^coeff(5))/(Final_Candidates{x_MolF,13}^coeff(6)),0);
                end
                Final_Candidates{x_MolF-No_Isomer,11}=No_Hist;    %Number of total histogram(s)
            end
            %%  Fragments in AP+ like [M]+
            Mass=Final_Matches(x(1),3);
            MolF=Chemical_structure_Print(ID_library(ID_i,:));
            if No_Isomer==1
                Final_Candidates{x_MolF,Column_xlx(1)}=MolF;
                Final_Candidates{x_MolF,Column_xlx(1)+1}=Mass;
            else
                Final_Candidates{x_MolF-No_Isomer,Column_xlx(1)}=MolF;
                Final_Candidates{x_MolF-No_Isomer,Column_xlx(1)+1}=Mass;
            end
            %% assigning candidate molecular formula based on molecular formula ions
            MolecularIon=ID_library(ID_i,:);
            c=MolecularIon(1);h=MolecularIon(2);br=MolecularIon(3);cl=MolecularIon(4);f=MolecularIon(5);i=MolecularIon(6);
            n=MolecularIon(7);na=MolecularIon(8);o=MolecularIon(9);p=MolecularIon(10);s=MolecularIon(11);
            %% Fragments in AP+ are like [M-Cl]+
            if cl>=1
                molecule=[c,h,br,cl+1,f,i,n,na,o,p,s];
                cl=molecule(4);
                x_XIC=find(c==XIC_Combination(:,1)...
                    & br==XIC_Combination(:,3)...
                    & cl==XIC_Combination(:,4)...
                    & s==XIC_Combination(:,11));
                if ~isempty(x_XIC)
                    B=XIC_Primary{x_XIC};
                    Most_abundant_mass_x= B(:,3)==100;
                    X=B(Most_abundant_mass_x,1);
                    Mass=round(Most_Abundant_Mass([c,h,br,cl,f,i,n,na,o,p,s],X),5);
                else
                    Mass='N/A';
                end
                MolF=Chemical_structure_Print(molecule);
                if No_Isomer==1
                    Final_Candidates{x_MolF,Column_xlx(2)}=MolF;
                    Final_Candidates{x_MolF,Column_xlx(2)+1}=Mass;
                else
                    Final_Candidates{x_MolF-No_Isomer,Column_xlx(2)}=MolF;
                    Final_Candidates{x_MolF-No_Isomer,Column_xlx(2)+1}=Mass;
                end
            end
            x_MolF=x_MolF+1;
            %%
            fprintf('%d - mass= %f \n\n',x_MolF-1,Final_Candidates{x_MolF-1,5})
        end
    end
end
Final_Candidates=Final_Candidates(1:x_MolF-1,:);
%% Checking the availabilty of the candidate compounds in the EPA DSS Tox database
load([route,'EPA_DSSTox.mat'])
Molocule_desalted=EPA_DSSTox(:,4);  % MS ready
Molocule_mixed=EPA_DSSTox(:,3);     % Might be in salt form or not
for i=1:size(Final_Candidates,1)
    if ~isempty(Final_Candidates{i,3})
        i
        for n_f=1:2         % n_f represent the number of ionization pathways
            if ~isempty(Final_Candidates{i,Column_xlx(n_f)})
                MolF=Final_Candidates{i,Column_xlx(n_f)};
                x_desalted=find(strcmp(Molocule_desalted,MolF)==1);
                x_mixed1=find(strcmp(Molocule_mixed,MolF)==1);
                x_mixed=setdiff(x_mixed1,x_desalted);
                No_Available_candidates=length(x_desalted)+length(x_mixed);
                Final_Candidates{i,Column_xlx(n_f)+2}=No_Available_candidates;
                if length(x_desalted)==1 && isempty(x_mixed)
                    Final_Candidates{i,Column_xlx(n_f)+3}=EPA_DSSTox{x_desalted,1};
                elseif length(x_mixed)==1 && isempty(x_desalted)
                    Final_Candidates{i,Column_xlx(n_f)+3}=EPA_DSSTox{x_mixed,1};
                end
            end
        end
    end
end
clear EPA_DSSTox Molocule_desalted Molocule_mixed
%% To check which compounds are already known
Final_Candidates=Known_Compounds_Search(Final_Candidates,Column_xlx);
%%
save('Final_Candidates.mat','Final_Candidates','-v7.3')
%% Making a spreadsheet from the results
k=1;
splitxls=700000;
while k<=size(Final_Candidates,1)
    namexls=['Final_Candidates_',num2str(k),'.xlsx'];
    xlswrite(namexls,{'ID in MATLAB code','Candidate molecular ion formula',...
        'Number of Isomer(s)','Detected mass of molecular ion (Da)','Exact mass of molecular ion',...
        'Intensity','Normalized Euclidean mass error (mDa)',['Profile Cosine Similarity (', char(8240),')'],...
        'Retention time (min)','Peak Area','NDCS','RCS(%)','RPW','Peak Identification score',...
        'Candidate molecular formula [M]+','Exact mass of candidate compound  [M]+ (Da)',...
        'EPA DSSTox','Name','Candidate molecular formula [M-Cl]+',...
        'Exact mass of candidate compound  [M-Cl]+ (Da)','EPA DSSTox','Name','Link'});
    if k+splitxls<=size(Final_Candidates,1)
        xlswrite(namexls,Final_Candidates(k:k+splitxls,:),['A2:V',num2str(splitxls+1)])
    else
        xlswrite(namexls,Final_Candidates(k:size(Final_Candidates,1),:),['A2:V',num2str(size(Final_Candidates,1)-k+1)])
    end
    k=k+splitxls
end
%% To create a folder containing of compounds that have at least one record in comptox
Chromatogram_Comptox_DSSTox_2
