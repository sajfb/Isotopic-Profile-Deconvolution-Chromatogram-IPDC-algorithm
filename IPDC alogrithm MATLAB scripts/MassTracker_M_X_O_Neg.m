%% This module is to compare candidate molecular formulas for a mass of interest
clc
mass_Target=484.70633;      % <-- Targted mass (Da)
RT_Target=44.818;           % <-- Targted retention time (min)
Delta_mass=0.01;            % Mass tolerance (Da)
Delta_RT=0.1;               % Retention time tolerance (min)
kmax=50;                    % Number of candidates
%% This section should be run once and when you have its variables in Workspace, you do not need to keep reruning this section
load('Final_Candidates.mat')
load('ID_library.mat')
EntireID_Isomer=Final_Candidates(:,1);
MassID_Isomer=Final_Candidates(:,4);
RTID_Isomer=Final_Candidates(:,9);
ScoreID_Isomer=Final_Candidates(:,14);
ID=Final_Candidates{1,1};
for i=1:size(EntireID_Isomer,1)
    if ~isnumeric(EntireID_Isomer{i})
        EntireID_Isomer{i}=ID;
    else
        ID=Final_Candidates{i,1};
        if Final_Candidates{i,3}~=1
            MassID_Isomer{i}=0;
            RTID_Isomer{i}=0;
            ScoreID_Isomer{i}=0;
        end
    end
end
EntireID_Isomer=cell2mat(EntireID_Isomer);
MassID_Isomer=cell2mat(MassID_Isomer);
RTID_Isomer=cell2mat(RTID_Isomer);
ScoreID_Isomer=cell2mat(ScoreID_Isomer);
%%
Counter=1;
A4=cell(kmax+1,21);
A4{Counter,3}=mass_Target;
A4{Counter,8}=RT_Target;
MplusOne=[0;1.00335;1.9974;3.0007;3.9948];%4.9981;5.9671;5.9922;6.9955;7.9896;8.9928];
x_Targeted1=[];
for M=1:length(MplusOne)
    x_Targeted1=[x_Targeted1;find(abs(mass_Target-MplusOne(M)-MassID_Isomer)<=Delta_mass & ...
        abs(RTID_Isomer-RT_Target)<Delta_RT)];
    x_Targeted1=[x_Targeted1;find(abs(mass_Target+MplusOne(M)-MassID_Isomer)<=Delta_mass & ...
        abs(RTID_Isomer-RT_Target)<Delta_RT)];
end
x_Targeted1=unique(x_Targeted1);
x_Targeted2=[x_Targeted1, ScoreID_Isomer(x_Targeted1),EntireID_Isomer(x_Targeted1)];
i=1;
while i<=size(x_Targeted2,1)
    x_f=find(x_Targeted2(i,2)==x_Targeted2(:,2) & x_Targeted2(i,3)==x_Targeted2(:,3));
    if length(x_f)>1
        x_Targeted2(i,:)=[];
        i=i-1;
    end
    i=i+1;
end
x_Targeted2=sortrows(x_Targeted2,-2);
x_Targeted=x_Targeted2(:,1);
A3=0;k=0;i=0;
while i<length(x_Targeted) && k<kmax
    i=i+1;
    A2=Final_Candidates(x_Targeted(i),4:14);
    if A2{8}>2 && A2{9}>25 && A2{10}<0.8 && k<kmax
        x_B=find(EntireID_Isomer==EntireID_Isomer(x_Targeted(i)));
        x_x_B=find(~cellfun(@isempty,Final_Candidates(x_B,19)));
        B=Final_Candidates(x_B(x_x_B),15:22);
        Counter=Counter+1;
        k=k+1;
        A1=['[',Chemical_structure_Print(ID_library(EntireID_Isomer(x_Targeted(i)),:)),']-'];
        A3=[k,A1,A2,B];
        A4(Counter,:)=A3;
    end
end
A4=A4(1:Counter,:);
%%
T = cell2table(A4,'VariableNames',...
    {'Order_of_candidate_compounds','Candidate_molecular_ion_formula','Detected_mass_of_molecular_ion',...
    'Exact_mass_of_molecular_ion','Intensity','Average_normalized_Euclidean_mass_error',...
    'Average_profile_cosine_similarity','Retention_time','Peak_Area','NDCS','RCS','RPW',...
    'Peak_Identification_score','Candidate_molecular_formula_M_X_O_',...
    'Exact_mass_of_candidate_compound_M_X_O_','Number_of_candidates_in_Comptox','Name',...
    'Candidate_molecular_formula_M_','Exact_mass_of_candidate_compound_M_','Number_of_candidates_in_Comptox2','Name2'})
