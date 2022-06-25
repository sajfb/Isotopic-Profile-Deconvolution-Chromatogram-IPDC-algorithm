delete(gcp)
clc;clearvars -except Peaks_TIC1
parpool('local',4)
% load('Peaks_TIC1.mat')
load('Retention_Time_TIC1.mat')
load('XIC_Quick.mat')
%% m/z screening input parameters
T0=6;                       % scan number of RT=6 min
T_end=55;                   % scan number of RT=55 min
Mass_error=15/1e3;          % mass error= 15mDa
% PPM1=20;                  % in order to use PPM for m/z screening section, lines 40 and 54 should be updated
min_INT=500;                % Intensity threshold
Min_Profile_CS=950;         % Profile cosine similarity
HIGHEST_mass=700;           % Maximum mass that the instrument scans or we are intersted in
%% Molecular formula assignment input parameters
Normalized_Mass_Error=10/1e3;   % Maximum allowed normalized EC mass error for entire isotopologues (Da)
% PPM2=20;                      % in order to use PPM, lines 150 and 162 should be updated
Mass_error_BasePeak=10/1e3;	% Maximum allowed mass error for the most abundant isotopologue (Da)
Min_Profile_CS2=950;            % Profile cosine similarity
%% m/z screening algorithm
[~,t0]=min(abs(Retention_Time_TIC1-T0));t0=t0(1);
[~,t_end]=min(abs(Retention_Time_TIC1-T_end));t_end=t_end(1);
S_XIC=size(XIC_Quick,1);
Primary_Match_dummy=cell(t_end-t0+1,2);
Counter2=0;
for t=t0:t_end
    tic
    t
    PEAKS=Peaks_TIC1{t};
    Q=zeros(S_XIC,1);
    G=cell(S_XIC,1);
    parfor i=1:S_XIC
        A_P=XIC_Quick{i};
        A_P1=A_P{2};
        MW_diff=A_P1(:,1);
        Theoretical_Intensity=A_P1(:,2);
        S=A_P{4};     % S=Number of isotopologues
        x_MAIso=A_P{5};    % MAIso(Most Abundant Isotopologue)
        k_initial=A_P{6};
        % Mass_error=PPM1*A_P(1,1)/1e6;
        m_z_x=find(PEAKS(:,1)>=A_P{1}-Mass_error & PEAKS(:,1)<HIGHEST_mass+8);
        M_Z_Exp=PEAKS(m_z_x,1);
        Intensity_Exp=PEAKS(m_z_x,2);
        D=[];
        Counter1=0;
        for j=1:length(M_Z_Exp)
            if Intensity_Exp(j)>min_INT
                index=cell(S,1);
                index{x_MAIso}=j;
                size_index=zeros(S,1);
                size_index(x_MAIso)=1;
                k=k_initial;
                while k<=S
                    % Mass_error=PPM1*M_Z_Exp(j)/1e6;
                    index{k,1}=find(abs(M_Z_Exp-M_Z_Exp(j)-MW_diff(k))<=Mass_error);
                    if isempty(index{k,1})
                        k=S+1;
                    else
                        size_index(k,1)=size(index{k},1);
                        if k+1~=x_MAIso
                            k=k+1;
                        else
                            k=x_MAIso+1;
                        end
                    end
                end
                P=prod(size_index);
                if P~=0
                    Candidate_Comparison=zeros(S,P);
                    for o=1:S
                        INDEX=index{o};
                        A=[];
                        k1=0;
                        while length(A)~=P
                            k1=k1+1;
                            A=[A,ones(1,P/prod(size_index(1:o)))*INDEX(k1)];
                            if k1==length(INDEX)
                                k1=0;
                            end
                        end
                        Candidate_Comparison(o,:)=A;
                    end
                    Comparison=cell(1,P);
                    Profile_CS=zeros(P,1);
                    B=zeros(S,2);
                    for n=1:P
                        B(:,1)=round(M_Z_Exp(Candidate_Comparison(:,n)),5);           % Experimental m/z
                        B(:,2)=round(Intensity_Exp(Candidate_Comparison(:,n)),0);     % Experimental Intensity
                        Comparison{1,n}=B;
                        Profile_CS(n)=sum(B(:,2).*Theoretical_Intensity)/sqrt(sum(B(:,2).^2))/sqrt(sum(Theoretical_Intensity.^2))*1000;   % Profile cosine similarity
                    end
                    [~,x_Primary]=max(Profile_CS);
                    if Profile_CS(x_Primary)>Min_Profile_CS
                        Counter1=Counter1+1;
                        D{Counter1,1}=Comparison{1,x_Primary};
                    end
                end
            end
        end
        if ~isempty(D)
            Q(i,1)=i;
            G{i,1}=D;
        end
    end
    x_Q=find(Q~=0);
    if ~isempty(x_Q)
        W=[num2cell(Q(x_Q)),G(~cellfun(@isempty,G))];
        Counter2=Counter2+1;
        Primary_Match_dummy{Counter2,1}=t;
        Primary_Match_dummy{Counter2,2}=W;
    end
    toc
end
clear Peaks_TIC1
delete(gcp)
Primary_Match=Primary_Match_dummy(1:Counter2,:);    % Primary_Match contains raw halogenated features 
% save('Primary_Match.mat','Primary_Match','-v7.3')
%%  Molecular formula assignment algorithm
clc;clearvars -except Primary_Match XIC_Quick Min_Profile_CS2 Normalized_Mass_Error Mass_error_BasePeak
% load('Primary_Match.mat')
load('XIC_Primary.mat')
load('ID_Mass.mat')
CounterF=0;
Final_Matches=zeros(10e7,7);
for t=1:size(Primary_Match,1)
    t
    ScanNumber=Primary_Match{t,1};
    T_PRIME_MATCH=Primary_Match{t,2};
    S_X=[];
    I_MATCH=[];
    Counter1=0;
    for i=1:size(T_PRIME_MATCH,1)
        S_X=[S_X;XIC_Quick{T_PRIME_MATCH{i,1}}{1,3}];
        for j=1:length(XIC_Quick{T_PRIME_MATCH{i,1}}{1,3})
            Counter1=Counter1+1;
            I_MATCH{Counter1,1}=T_PRIME_MATCH{i,2};
        end
    end
    for i=1:size(S_X,1)
        II_MATCH=I_MATCH{i};
        S_MATCH=size(II_MATCH,1);
        B=XIC_Primary{S_X(i),1};
        S=XIC_Primary{S_X(i),2};
        x_100=XIC_Primary{S_X(i),3};
        ID_initial=XIC_Primary{S_X(i),4};
        ID_end=XIC_Primary{S_X(i),5};
        for kop=1:S_MATCH
            A=II_MATCH{kop};
            m_z=A(x_100,1);
            % Mass_error_BasePeak=PPM2*Mass/1e6;
            C1=find(abs(ID_Mass(ID_initial:ID_end,1)-m_z)<=Mass_error_BasePeak);
            if ~isempty(C1)
                C1=ID_initial-1+C1;
                A(:,3)=A(:,2)/A(x_100,2)*100;
                A(:,4)=B(:,3);
                for iC1=1:size(C1,1)
                    HOFNPNa_mass=ID_Mass(C1(iC1),2);
                    Exact_mass=round(B(:,1)+HOFNPNa_mass,5);
                    mass_error=abs(Exact_mass-A(:,1));
                    Profile_CS=sum(A(:,3).*A(:,4))/sqrt(sum(A(:,3).^2))/sqrt(sum(A(:,4).^2))*1000;   % Profile cosine similarity
                    EUC_Msss_Error=sqrt(sum(mass_error.^2)/S);             % Normalized euclidean distance for mass error (mDa)
                    % EUC_Msss_Error=PPM2*Exact_mass(1,1)/1e6;
                    if EUC_Msss_Error<=Normalized_Mass_Error && Profile_CS>Min_Profile_CS2
                        CounterF=CounterF+1;
                        Final_Matches(CounterF,1)=C1(iC1);
                        Final_Matches(CounterF,2)=round(A(x_100,1),5);                      % Experimental Mass
                        Final_Matches(CounterF,3)=round(Exact_mass(x_100,1),5);             % Exact Mass
                        Final_Matches(CounterF,4)=A(x_100,2);                               % Intensity
                        Final_Matches(CounterF,5)=round(EUC_Msss_Error*1e3,2);              % in mDa
                        Final_Matches(CounterF,6)=round(Profile_CS,3);                      % in per-mille
                        Final_Matches(CounterF,7)=ScanNumber;
                    end
                end
            end
        end
    end
end
clc;clearvars -except Final_Matches CounterF
Final_Matches=Final_Matches(1:CounterF,:);
Final_Matches=sortrows(Final_Matches,7);   % To sort the matches based on the their scan number
Final_Matches=sortrows(Final_Matches,1);   % To sort the matches based on the their ID
Final_Matches=sortrows(Final_Matches,-3);  % To sort the matches based on the their masses
save('Final_Matches.mat','Final_Matches','-v7.3')
%%
