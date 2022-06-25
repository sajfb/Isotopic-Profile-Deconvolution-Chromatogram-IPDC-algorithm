function [PA_total, No_Isomer]=PeakAreaPlot(Peaks_TIC1,Retention_Time_TIC1,MZ,Mass_window,Int,Scan_Number,...
    Chromatogram_scan_interval,Smoothing_window,noise_removal,Max_No_isomers)
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
hold on
box on
Missed_peaks=0;
Enitre_SN=[];
Segment=Peak_detection(Chromatogram,Chromatogram_scan_interval,noise_removal,Max_No_isomers);
No_Isomer=size(Segment,1);
ScanNumber_isomers=cell(No_Isomer,1);
PA_isomers=zeros(No_Isomer,1);
for j1=1:No_Isomer
    if Segment(j1,1)==1
        Segment(j1,1)=Segment(j1,1)+1;
    end
    Seg=Segment(j1,1):Segment(j1,2);
    ScanNumber_isomers{j1}=intersect(Scan_Number,Seg);
    Enitre_SN=[Enitre_SN;ScanNumber_isomers{j1}];
    if ~isempty(ScanNumber_isomers{j1})
        Seg1=find(Chromatogram(:,1)==Segment(j1,1)); % beginning of the isomer peak
        Seg2=find(Chromatogram(:,1)==Segment(j1,2)); % end of the isomer peak
        chromW=Chromatogram(Seg1:Seg2,:);
        [~,Seg1_RAW]=min(abs(ChromatogramRAW(:,1)-Chromatogram(Seg1,1)));
        [~,Seg2_RAW]=min(abs(ChromatogramRAW(:,1)-Chromatogram(Seg2,1)));
        max_C_W=max(chromW(:,2));
        max_C_RAW=max(ChromatogramRAW(Seg1_RAW:Seg2_RAW,2));
        chromW(:,2)=chromW(:,2)*max_C_RAW/max_C_W;
        if Chromatogram(Seg1,2)~=0
            chromW=[[chromW(1,1)-1,0];chromW];
        end
        if Chromatogram(Seg2,2)~=0
            chromW=[chromW;[chromW(size(chromW,1),1)+1,0]];
        end
        PA_isomers(j1)=trapz(Retention_Time_TIC1(chromW(:,1)),chromW(:,2));
        plot(Retention_Time_TIC1(chromW(:,1)),chromW(:,2),'LineStyle',':','LineWidth',2.5)
    else
        Missed_peaks=Missed_peaks+1;
    end
end
if Missed_peaks>0
    PA_isomers=PA_isomers(PA_isomers~=0);
    No_Isomer=No_Isomer-Missed_peaks;
end
PA_total=round(sum(PA_isomers),2);
x_INT1=intersect(Enitre_SN,Scan_Number);
x_INT=zeros(size(Scan_Number,1),1);
for i=1:size(Scan_Number,1)
    if ismember(Scan_Number(i),x_INT1)
        x_INT(i,1)=i;
    end
end
x_INT=x_INT(x_INT~=0);
x_INT_RT=Scan_Number(x_INT);
stem(Retention_Time_TIC1(x_INT_RT),Int(x_INT),'Marker','none','Color',[0 0.447058826684952 0.74117648601532])
T=[Retention_Time_TIC1(x_INT_RT(1))-0.25, Retention_Time_TIC1(x_INT_RT((size(x_INT_RT,1))))+0.25];
xlim(T)
sinle_chrom_mz=CallmzPeaks_TIC1(MZ,T,Mass_window,Peaks_TIC1,Retention_Time_TIC1);
plot(sinle_chrom_mz(:,2),sinle_chrom_mz(:,1),'LineWidth',0.5,'LineStyle','-.','Color',[0 0 0]);
strPA=['Peak Area = ',num2str(round(PA_total,0))];
annotation('textbox',[0.747356051703878 0.954057170577345 0.158112698296122 0.0349492664872055],...
    'String',strPA,...
    'LineStyle','none',...
    'FontSize',16);
xlabel('Retention Time (min)');
ylabel('Intensity');
hold off
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
function sinle_chrom_mz=CallmzPeaks_TIC1(mz,T,Mass_window,Peaks_TIC1,Retention_Time_TIC1)
T0=T(1);Tend=T(2);
[~,t0]=min(abs(Retention_Time_TIC1-T0));t0=t0(1);
[~,t_end]=min(abs(Retention_Time_TIC1-Tend));t_end=t_end(1);
sinle_chrom_mz=zeros(t_end-t0+1,2);
for t=t0:t_end
    PEAKS=Peaks_TIC1{t};
    x=find(abs(PEAKS(:,1)-mz)<=Mass_window);
    if ~isempty(x)
        if length(x)>1
            [~,y]=min(abs(PEAKS(x,1)-mz));
            x=x(y(1));
        end
        sinle_chrom_mz(t,:)=[PEAKS(x,2),Retention_Time_TIC1(t)];
    else
        sinle_chrom_mz(t,2)=Retention_Time_TIC1(t);
    end
end
