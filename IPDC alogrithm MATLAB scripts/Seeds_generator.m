clear;clc;
%%
LOWEST_mass=50;         % Minimum mass that the instrument scans or we are intersted in
HIGHEST_mass=700;       % Maximum mass that the instrument scans or we are intersted in
%%
Peak_Spacing=1e-3;      % Based on the instrument resolution to resolve neighboring isotopologues %(Da)
% The_Int_Cutoff=6;      % Theoretical intensity threshold to cut off low abundant isotoplogues
                         % In this work, this value is set according to number of carbons in line 22
%% To combine similar seed profiles
MAMD=1e-3;      % Maximum Allowed Mass differences %(Da)
MPCS=999;       % Minimum of Profile Cosine Similarity (per-mille)
%%
counter=0;
SIZE=1; % Expected size of the isotope pattern seeds (OPTIONAL VALUE)
XIC_Combination=zeros(SIZE,11);
XIC_Primary=cell(SIZE,5);
N=zeros(SIZE,1);
ID_library=zeros(10e7,11);
ID_Mass=zeros(10e7,2);
Counter_ID=0;
for c=4:25  % Carbon range
    The_Int_Cutoff=min(10,c);
    for br=0:8  % Bromine range
        for cl=0:12  % Chlorine range
            if br+cl>=1 && br+cl<=16    % Optional condition to remove low probable compounds
                for s=0:2  % Sulfur range
                    Molecule=[c 0 br cl 0 0 0 0 0 0 s];
                    [MW_iso, intensity_iso]=Isotopic_Profile(Molecule,Peak_Spacing,The_Int_Cutoff);
                    if max(MW_iso)<=HIGHEST_mass
                        counter=counter+1
                        XIC_Combination(counter,:)=Molecule;
                        x_100=find(intensity_iso==100);
                        MW_diff=MW_iso-MW_iso(x_100);
                        N(counter,1)=size(MW_diff,1);
                        XIC_Primary{counter,1}=[round(MW_iso,5), MW_diff, round(intensity_iso,2)];
                        XIC_Primary{counter,2}=N(counter,1);
                        XIC_Primary{counter,3}=x_100;
                        %% ID library
                        XIC_Primary{counter,4}=Counter_ID+1;
                        for h=0:2*c+1  % Hydrogen range
                            for f=0 % Flourine range
                                for i=0:0   % Iodine range
                                    for n=0:3  % Nitrogen range
                                        if h+cl+br+f+i>=c/2-2*n % Optional condition to remove low probable compounds (rule 5)
                                            for na=0:0  % Sodium range
                                                for o=0:5  % Oxygen range
                                                    for p=0:1  % Phosphorus range
                                                        MASS=Most_Abundant_Mass([c,h,br,cl,f,i,n,na,o,p,s],MW_iso(x_100));
                                                        if MASS<=HIGHEST_mass && MASS>=LOWEST_mass
                                                            Counter_ID=Counter_ID+1;
                                                            ID_library(Counter_ID,:)=[c,h,br,cl,f,i,n,na,o,p,s];
                                                            ID_Mass(Counter_ID,1)=MASS;
                                                            ID_Mass(Counter_ID,2)=MASS-MW_iso(x_100);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        XIC_Primary{counter,5}=Counter_ID;
                        %%
                    end
                end
            end
        end
    end
end
ID_library=ID_library(1:Counter_ID,:);
ID_Mass=ID_Mass(1:Counter_ID,:);
save('XIC_Combination.mat','XIC_Combination')
save('XIC_Primary.mat','XIC_Primary')
save('ID_Mass.mat','ID_Mass','-v7.3')
save('ID_library.mat','ID_library','-v7.3')
% % % % % % % % xlswrite('XIC_Combination.xlsx',XIC_Combination)
%% To combine similar seed profiles
counterQ=0;
SIZE=1; % Expected size of the contracted isotope pattern seeds (OPTIONAL VALUE)
XIC_Quick=cell(SIZE,1);
for i=1:size(XIC_Primary,1)
    i
    ID_P=[];   % ID_P is to record the number ID of the similar patterns
    if N(i)~=0
        S=N(i);  % number of isotopologues
        n=find(S==N);
        A=XIC_Primary{i,1};
        min_Mass=A(1,1);
        counter=0;
        for kn=n'
            B=XIC_Primary{kn,1};
            C=abs(A-B);
            Max_M_D=max(C(:,2));    % Mass difference
            PCS=sum(A(:,3).*B(:,3))/sqrt(sum(A(:,3).^2))/sqrt(sum(B(:,3).^2))*1000;  % Profile cosine similarity
            if Max_M_D<=MAMD && PCS>=MPCS
                counter=counter+1;
                if B(1,1)<min_Mass
                    min_Mass=B(1,1);
                end
                ID_P(counter,:)=kn;
                N(kn)=0;
            end
        end
    end
    if ~isempty(ID_P)
        Entire_mass_diff=zeros(S,length(ID_P));
        Entire_profile=zeros(S,length(ID_P));
        for p=1:length(ID_P)
            TTT=XIC_Primary{ID_P(p,1),1};
            Entire_mass_diff(:,p)=TTT(:,2);
            Entire_profile(:,p)=TTT(:,3);
        end
        Ave_mass_diff=zeros(S,1);
        Ave_profile=zeros(S,1);
        for k=1:S
            Ave_mass_diff(k)=round(mean(Entire_mass_diff(k,:)),5);
            Ave_profile(k)=round(mean(Entire_profile(k,:)),2);
        end
        Profile_pattern=[Ave_mass_diff Ave_profile];
        counterQ=counterQ+1;
        S=size(Profile_pattern,1);     % S=Number of isotopologues
        [~, x_MAIso]=max(Profile_pattern(:,2));    % MAIso(Most Abundant Isotopologue)
        if x_MAIso~=1    % if the most abundant isotopologue is the monoisotopic isotopologues
            k_initial=1;
        else
            k_initial=2;
        end
        XIC_Quick{counterQ,1}={min_Mass, Profile_pattern, ID_P, S, x_MAIso, k_initial};
    end
end
save('XIC_Quick.mat','XIC_Quick')
