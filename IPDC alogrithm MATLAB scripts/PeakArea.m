function [No_Isomer, ScanNumber_isomers, PA_isomers, n_w_Int, RPW, neme,...
    pcs]=PeakArea(Peaks_TIC1,Retention_Time_TIC1,IsotopeModel,Mass_error,...
    Int,Scan_Number,Chromatogram_scan_interval,Smoothing_window,...
    noise_removal,Max_No_isomers,min_INT)
Chromatogram=zeros(Scan_Number(size(Int,1))-Scan_Number(1)+1,2);
for a=1:Scan_Number(size(Int,1))-Scan_Number(1)+1
    index1=find(Scan_Number==Scan_Number(1)+a-1, 1);
    if ~isempty(index1)
        Chromatogram(a,:)=[Scan_Number(1)+a-1 Int(index1(1))];
    else
        Chromatogram(a,:)=[Scan_Number(1)+a-1 0];
    end
end
ChromatogramRAW=Chromatogram;
%% Filling missing chromatogram scans
SZC=size(Chromatogram,1);
i=1;
while i<=SZC
    if Chromatogram(i,2)~=0
        if i+Chromatogram_scan_interval-1>=SZC
            i_end=SZC;
        else
            i_end=i+Chromatogram_scan_interval-1;
        end
        x=find(Chromatogram(i:i_end,2)>0);
        i_f=i;
        while ~isempty(x)
            if i+Chromatogram_scan_interval-1>=SZC
                i_end=SZC;
            else
                i_end=i+Chromatogram_scan_interval-1;
            end
            x=find(Chromatogram(i:i_end,2)>0);
            if ~isempty(x)
                i=i+x(length(x));
            else
                i=i-1;
                i_e=i;
                i_e_f=i_f:i_e;
                if length(i_e_f)~=1
                    index=find(Chromatogram(i_e_f,2)~=0);
                    x_fill=Chromatogram(i_e_f(index),1);
                    v_fill=Chromatogram(i_e_f(index),2);
                    xq_fill=Chromatogram(i_f,1):Chromatogram(i_e,1);
                    Chromatogram(i_f:i_e,2)=interp1(x_fill,v_fill,xq_fill,'phcip');
                end
                break
            end
        end
    end
    i=i+1;
end
%% Smoothing the chromatogram trace
C1=[(Chromatogram(1,1)-Chromatogram_scan_interval:Chromatogram(1,1)-1)',zeros(Chromatogram_scan_interval,1)];
C1(C1(:,1)<=1,:)=[];
C2=[(Chromatogram(SZC,1)+1:Chromatogram(SZC,1)+Chromatogram_scan_interval)',zeros(Chromatogram_scan_interval,1)];
Chromatogram=[C1;Chromatogram;C2];
Chromatogram(:,2)=smoothdata(Chromatogram(:,2),'lowess',Smoothing_window,'omitnan');
Chromatogram(Chromatogram(:,2)<0,2)=0;
Chromatogram(Chromatogram(:,1)>=length(Retention_Time_TIC1),:)=[];
%%
Missed_peaks=0;
Segment=Peak_detection(Chromatogram,Chromatogram_scan_interval,noise_removal,Max_No_isomers);
No_Isomer=size(Segment,1);
Chromatogram=[[Chromatogram(1,1)-1,0];Chromatogram;[Chromatogram(size(Chromatogram,1),1)+1,0]];
ScanNumber_isomers=cell(No_Isomer,1);
PA_isomers=zeros(No_Isomer,1);
RPW=zeros(No_Isomer,1);   % Percentage of Guassianity Height 
n_w_Int=zeros(No_Isomer,1);
for j1=1:No_Isomer
    if Segment(j1,1)==1
        Segment(j1,1)=Segment(j1,1)+1;
    end
    Seg=Segment(j1,1):Segment(j1,2);
    ScanNumber_isomers{j1}=intersect(Scan_Number,Seg);
    if ~isempty(ScanNumber_isomers{j1})
        Seg1=find(Chromatogram(:,1)==Segment(j1,1)); % beginning of the isomer peak
        Seg2=find(Chromatogram(:,1)==Segment(j1,2)); % end of the isomer peak
        chromW=Chromatogram(Seg1:Seg2,:);
        %%
        [~,Seg1_RAW]=min(abs(ChromatogramRAW(:,1)-Chromatogram(Seg1,1)));
        [~,Seg2_RAW]=min(abs(ChromatogramRAW(:,1)-Chromatogram(Seg2,1)));
        max_C_W=max(chromW(:,2));
        max_C_RAW=max(ChromatogramRAW(Seg1_RAW:Seg2_RAW,2));
        chromW(:,2)=chromW(:,2)*max_C_RAW/max_C_W;
        %%
        if Chromatogram(Seg1,2)~=0
            chromW=[[chromW(1,1)-1,0];chromW];
        end
        if Chromatogram(Seg2,2)~=0
            chromW=[chromW;[chromW(size(chromW,1),1)+1,0]];
        end
        [neme(j1),pcs(j1)]=CombineSpectra(Peaks_TIC1,IsotopeModel,chromW(1,1),chromW(size(chromW,1),1),Mass_error);
        PA_isomers(j1)=round(trapz(Retention_Time_TIC1(chromW(:,1)),chromW(:,2)),2);
        [n_w_Int(j1), RPW(j1)]=PeakWidthCalculator(Retention_Time_TIC1,chromW,min_INT);
    else
        Missed_peaks=Missed_peaks+1;
    end
end
if Missed_peaks>0
    ScanNumber_isomers=ScanNumber_isomers(~cellfun(@isempty,ScanNumber_isomers));
    PA_isomers=PA_isomers(PA_isomers~=0);
    n_w_Int=n_w_Int(n_w_Int~=0);
    RPW=RPW(RPW~=0);
    No_Isomer=No_Isomer-Missed_peaks;
end
%%
function Segments=Peak_detection(Chromatogram,Chromatogram_scan_interval,noise_removal,Max_No_isomers)
Scan_Number=Chromatogram(:,1);
B=[Chromatogram(:,1),Chromatogram(:,2),-double(islocalmin(Chromatogram(:,2)))];
% To merge close local minimums
SZC=size(Chromatogram(:,1),1);
x=[];
for i=1:SZC
    if B(i,3)==-1 && i<=Chromatogram_scan_interval
        x=find(B(1:i,3)==-1);
    elseif B(i,3)==-1 && SZC-i<=Chromatogram_scan_interval
        x=find(B(i:SZC,3)==-1);
        x=x+i-1;
    elseif B(i,3)==-1
        x=find(B(i-Chromatogram_scan_interval:i+Chromatogram_scan_interval,3)==-1);
        x=x+i-Chromatogram_scan_interval-1;
    end
    if ~isempty(x) && size(x,1)>1
        y=find(min(B(x,2))==B(x,2) & B(x,2)~=0);
        if ~isempty(y)
            B(x,3)=0;
            B(x(y(1)),3)=-1;
        end
        x=[];
    end
end
x=find(B(:,3)==-1);
Int_0=find(B(:,2)==0);
Segments=[[1;x],[x;size(B,1)]];
for i=1:size(Segments,1)
    x=find(Segments(i,1)<=Int_0 & Segments(i,2)>=Int_0);
    if ~isempty(x)
        L=setdiff(Segments(i,1):Segments(i,2),Int_0(x));
        if ~ismember(Segments(i,1),L)
            Segments(i,1)=L(1)-1;
        end
        if ~ismember(Segments(i,2),L)
            Segments(i,2)=max(L)+1;
        end
    end
end
%% To merge the fornt or tail of the peaks
x_S1=Segments(:,1);
x_S2=Segments(:,2);
x_max_Seg=zeros(size(Segments,1),1);
for i=1:size(Segments,1)
    x_max_Seg(i,1)=find(max(Chromatogram(x_S1(i):x_S2(i),2))==Chromatogram(:,2),1);
end
R=[x_S1 x_max_Seg x_S2];
i=1;
while i~=size(R,1)
    if R(i+1,1)-R(i,3)<=Chromatogram_scan_interval
        r2=Chromatogram(R(i,2),1);
        r3=Chromatogram(R(i,3),1);
        r5=Chromatogram(R(i+1,2),1);
        l2=Chromatogram(R(i,2),2);
        l5=Chromatogram(R(i+1,2),2);
        if l2>=l5
            l3=l5-Chromatogram(R(i+1,1),2);
            z=(l2-l5)*(r5-r3)/(r5-r2);
        elseif l2<l5
            l3=l2-Chromatogram(R(i,1),2);
            z=(l5-l2)*(r3-r2)/(r5-r2);
        end
        if l3/z<5e-2
            if l2>=l5
                R(i,3)=R(i+1,3);
                R(i+1,:)=[];
                i=i-1;
            elseif l2<l5
                R(i,:)=[R(i,1) R(i+1,2) R(i+1,3)];
                if i~=size(R,1)
                    R(i+1,:)=[];
                    i=i-1;
                end
            end
        end
    end
    i=i+1;
end
Segments=[Scan_Number(R(:,1)),Scan_Number(R(:,3))];
%%
S_N1=zeros(size(Segments,1),1);
for i=1:size(Segments,1)
    x_S_N1=find(Chromatogram(:,1)==Segments(i,1));
    x_S_N2=find(Chromatogram(:,1)==Segments(i,2));
    S_N1(i)=max(Chromatogram(x_S_N1:x_S_N2,2));
end
x_S_N1=find(S_N1/max(S_N1)*100>noise_removal);  % noise level below %
if size(x_S_N1,1)>=Max_No_isomers
    [~,x_S_N_k]=maxk(S_N1,Max_No_isomers);
    x_S_N1=intersect(x_S_N_k,x_S_N1);
end
Segments=Segments(x_S_N1,:);
%%
function [n_w_Int, RPW]=PeakWidthCalculator(Retention_Time_TIC1,chromW,min_INT)
[H, x_H]=max(chromW(:,2));
S=size(chromW,1);
[V_mode, f]=mode(chromW(1:x_H,2));
while f~=1
    Ind=1:x_H;
    x=find(chromW(Ind,2)==V_mode);
    for k=2:size(x,1)
        chromW(Ind(x(k)),2)=chromW(Ind(x(k)),2)+0.0001*k;
    end
    [V_mode, f]=mode(chromW(1:x_H,2));
end
[V_mode, f]=mode(chromW(x_H:S,2));
while f~=1
    Ind=x_H:S;
    x=find(chromW(Ind,2)==V_mode);
    for k=2:size(x,1)
        chromW(Ind(x(k)),2)=chromW(Ind(x(k)),2)+0.0001*k;
    end
    [V_mode, f]=mode(chromW(x_H:S,2));
end
Peak_Width=Retention_Time_TIC1(chromW(S,1))-Retention_Time_TIC1(chromW(1,1));
%%
RT_500_F=interp1(chromW(1:x_H,2),Retention_Time_TIC1(chromW(1:x_H,1)),min_INT);
RT_500_E=interp1(chromW(x_H:S,2),Retention_Time_TIC1(chromW(x_H:S,1)),min_INT);
x_RT_500_F=find((Retention_Time_TIC1-RT_500_F)>0,1);
x_RT_500_E=find((Retention_Time_TIC1-RT_500_E)<0);x_RT_500_E=x_RT_500_E(length(x_RT_500_E));
n_w_Int=x_RT_500_E-x_RT_500_F+1;   % Number of scans at 500
%%
RT_H_1_2_F=interp1(chromW(1:x_H,2),Retention_Time_TIC1(chromW(1:x_H,1)),H/2);    % index of H1/2
RT_H_1_2_E=interp1(chromW(x_H:S,2),Retention_Time_TIC1(chromW(x_H:S,1)),H/2);
Peak_Width_1_2=RT_H_1_2_E-RT_H_1_2_F;   % half peak width
RPW=Peak_Width_1_2/Peak_Width;
function [neme,pcs]=CombineSpectra(Peaks_TIC1,IsotopeModel,t0,t_end,MassError)
Theoretical_MW=IsotopeModel(:,1);
Theoretical_Intensity=IsotopeModel(:,2);
S=size(Theoretical_MW,1);
A=cell(t_end-t_end+1,S);
for t=t0:t_end
    PEAKS=Peaks_TIC1{t};
    M_Z_Exp=PEAKS(:,1);
    for i=1:S
        x=find(abs(M_Z_Exp-Theoretical_MW(i))<MassError);
        if length(x)>1
            [~,x_min_mass]=min(abs(M_Z_Exp-Theoretical_MW(i)));
            x=x_min_mass(1);
        end
        A{t-t0+1,i}=PEAKS(x,:);
    end
end
m_z_cluster=zeros(S,1);Int_cluster=zeros(S,1);
for i=1:S
    WS=cell2mat(A(:,i));
    Int_cluster(i,1)=sum(WS(:,2));
    m_z_cluster(i,1)=sum(WS(:,1).*WS(:,2))/Int_cluster(i);
end
neme=sqrt(sum((m_z_cluster-Theoretical_MW).^2)/S)*1000;
pcs=sum(Int_cluster.*Theoretical_Intensity)/sqrt(sum(Int_cluster.^2)*sum(Theoretical_Intensity.^2))*1000;   % Profile cosine similarity
