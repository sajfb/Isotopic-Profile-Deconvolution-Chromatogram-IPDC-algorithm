%% Molecular_Formula=[Carbon Hydrogen Bromine Chlorine Fluorine Idoine Nitrogen Sodium Oxygen Phosphorus Sulfur]
% Satellite group model is applied instead of rectangular model to merge neighboring isotopologues
% Peak_Spacing=1e-3;    % To merge neighboring isotopologues that can not be resolved within the instrument resolution(mDa)
% Molecular_Formula=[8,0,0,0,17,0,0,0,3,0,1];   % molecular formula vector for [C8F17O3S]- PFOS
% intensity_cutoff=4;                           % Intensity cutoff for relative abundancy (%)
function [MW_filtered, intensity_filtered]=Isotopic_Profile(Molecular_Formula,Peak_Spacing,intensity_cutoff)
C=[12 13.00335];
RA_C=[98.93 1.07]/100; %Relative abundancy of C
H=[1.00783 2.01410];
RA_H=[99.9885 0.0115]/100; %Relative abundancy of H
Br=[78.91834 80.91629];
RA_Br=[50.69 49.31]/100; %Relative abundancy of Br
Cl=[34.96885 36.96590];
RA_Cl=[75.78 24.22]/100; %Relative abundancy of Cl
F=18.9984032;
RA_F=100/100; %Relative abundancy of F
I=126.90448;
RA_I=100/100; %Relative abundancy of I
N=[14.00307 15.00011];
RA_N=[99.632 0.368]/100; %Relative abundancy of N
Na=22.98976928;
RA_Na=100/100; %Relative abundancy of Na
O=[15.99491 16.99913 17.99916];
RA_O=[99.757 0.038 0.205]/100; %Relative abundancy of O
P=30.97377;
RA_P=100/100; %Relative abundancy of P
S=[31.97207 32.97146 33.96885 35.96708];
RA_S=[94.93 0.76 4.29 0.02]/100; %Relative abundancy of S
%%
[C_combination,RAc_C]=SUB_COMB(Molecular_Formula(1),C,RA_C);
[H_combination,RAh_H]=SUB_COMB(Molecular_Formula(2),H,RA_H);
[Br_combination,RAbr_Br]=SUB_COMB(Molecular_Formula(3),Br,RA_Br);
[Cl_combination,RAcl_Cl]=SUB_COMB(Molecular_Formula(4),Cl,RA_Cl);
[F_combination,RAf_F]=SUB_COMB(Molecular_Formula(5),F,RA_F);
[I_combination,RAi_I]=SUB_COMB(Molecular_Formula(6),I,RA_I);
[N_combination,RAn_N]=SUB_COMB(Molecular_Formula(7),N,RA_N);
[Na_combination,RAn_Na]=SUB_COMB(Molecular_Formula(8),Na,RA_Na);
[O_combination,RAo_O]=SUB_COMB(Molecular_Formula(9),O,RA_O);
[P_combination,RAn_P]=SUB_COMB(Molecular_Formula(10),P,RA_P);
[S_combination,RAs_S]=SUB_COMB(Molecular_Formula(11),S,RA_S);
%%
Combination_Size=size(C_combination,1)*size(H_combination,1)*size(Br_combination,1)*...
    size(Cl_combination,1)*size(F_combination,1)*size(I_combination,1)*size(N_combination,1)*...
    size(Na_combination,1)*size(O_combination,1)*size(P_combination,1)*size(S_combination,1);
MW=zeros(Combination_Size,1);
RA=ones(Combination_Size,1); %RA indicates the abundance of isotopic combinations
Counter1=0;
for c=1:size(C_combination,1)  %number of C_combination
    for h=1:size(H_combination,1) %number of H_combination
        for br=1:size(Br_combination,1) %number of Br_combination
            for cl=1:size(Cl_combination,1) %number of Cl_combination
                for f=1:size(F_combination,1) %number of F_combination
                    for i=1:size(I_combination,1) %number of I_combination
                        for n=1:size(N_combination,1) %number of N_combination
                            for na=1:size(Na_combination,1) %number of Na_combination
                                for o=1:size(O_combination,1) %number of O_combination
                                    for p=1:size(P_combination,1) %number of P_combination
                                        for s=1:size(S_combination,1) %number of S_combination
                                            Counter1=Counter1+1;
                                            MW(Counter1)=sum(C_combination(c,:))+sum(H_combination(h,:))+...
                                                sum(Br_combination(br,:))+sum(Cl_combination(cl,:))+...
                                                sum(F_combination(f,:))+sum(I_combination(i,:))+...
                                                sum(N_combination(n,:))+sum(Na_combination(na,:)+...
                                                sum(O_combination(o,:))+sum(P_combination(p,:)+sum(S_combination(s,:))));
                                            RA(Counter1)=RAc_C(c)*RAh_H(h)*RAbr_Br(br)*RAcl_Cl(cl)*RAf_F(f)*...
                                                RAi_I(i)*RAn_N(n)*RAn_Na(na)*RAo_O(o)*RAn_P(p)*RAs_S(s);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
RA=RA*100;
B=[MW RA];
B=B(B(:,2)>1e-5,:);
B=sortrows(B,-2);
Counter=0;
if Peak_Spacing~=0
    for i=1:size(B,1)
        if B(i,1)~=0
            x=find(B(:,1)<B(i,1)+Peak_Spacing & B(:,1)>B(i,1)-Peak_Spacing);
            Counter=Counter+1;
            RA2(Counter,1)=sum(B(x,2));
            MW2(Counter,1)=sum(B(x,1).*B(x,2))/RA2(Counter,1);
            B(x,:)=0;
        end
    end
else
    RA2(:,1)=B(:,2);
    MW2(:,1)=B(:,1);
end
intensity=RA2/max(RA2)*100; %intensity or relative abundance
A=[MW2, intensity];
A=sortrows(A,1);
if Molecular_Formula(11)>0
    if Molecular_Formula(3)==0 && Molecular_Formula(4)==0
        intensity_cutoff=4;
    end
end
INDEX=find(A(:,2)>=intensity_cutoff);
if size(INDEX,1)==2
    M=sortrows(A(:,2),-1);
    m=find(M(2)==A(:,2));
    INDEX=[1;m(1)];
end
A=A(INDEX,:);
MW_filtered=A(:,1);intensity_filtered=A(:,2);
end
function combs = nmultichoosek(values, k)
%% // Return number of multisubsets or actual multisubsets.
%% http://stackoverflow.com/questions/28284671/generating-all-combinations-with-repetition-using-matlab
if numel(values)==1
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
end
function [element_combination,RAe_element]=SUB_COMB(Number_atoms,Isotope_masses,Isotope_Compositions)
%%// Subfunction for element and isotopic combinations.
if Number_atoms==0
    element_combination=0;
    RAe_element=1;
else
    if length(Isotope_Compositions)==1
        element_combination=Isotope_masses*ones(1,Number_atoms);
        RAe_element=Isotope_Compositions*ones(1,Number_atoms);
    else
        element_combination=nmultichoosek(Isotope_masses,Number_atoms);
        RAe_element=SUB_ISO_PRO(Number_atoms,Isotope_masses,Isotope_Compositions);
        x_c=find(RAe_element>1e-20);
        RAe_element=RAe_element(x_c);
        element_combination=element_combination(x_c,:);
    end
end
end
function Ai=SUB_ISO_PRO(number_of_Molecular_Formulas,Molecular_Formulaic_Mass,Isotopic_Composition)
Molecular_Formulaic_Mass_Combination=nmultichoosek(Molecular_Formulaic_Mass,number_of_Molecular_Formulas);
Ai=zeros(size(Molecular_Formulaic_Mass_Combination,1),1);
for k=1:size(Molecular_Formulaic_Mass_Combination,1)
    Combination=Molecular_Formulaic_Mass_Combination(k,:);
    F=1;
    R=1;
    for i=1:size(Molecular_Formulaic_Mass,2)
        x=find(Combination==Molecular_Formulaic_Mass(i));
        X=size(x,2);
        F=F*factorial(X);
        R=R*(Isotopic_Composition(i)^X);
    end
    Ai(k,1)=factorial(number_of_Molecular_Formulas)/F*R;
end
end
