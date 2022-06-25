%% This algorithm is designed to convert MS^e structure files into the MATLAB readable format
clear;clc;
mzxml_struct = mzxmlread('?????.mzxml');%==>Please input the filename here
Retention_Time_str={mzxml_struct.scan.retentionTime}';
basePeakIntensity=[mzxml_struct.scan.basePeakIntensity]';
basePeakMz=[mzxml_struct.scan.basePeakMz]';
ScanType={mzxml_struct.scan.scanType}';
msLEVEL=[mzxml_struct.scan.msLevel]';
size_Spectra=size(ScanType,1);
Retention_Time=zeros(size_Spectra,1);
peaks=cell(size_Spectra,1);
for i=1:size_Spectra
    Retention_Time(i)=sscanf(Retention_Time_str{i},'PT%f');
    A=mzxml_struct.scan(i).peaks.mz;
    i_m_z=1:2:length(A)-1;
    i_int=2:2:length(A);
    peaks{i}=[A(i_m_z), A(i_int)];
end
%% Importing msLevel 1 peaks
index0=find(basePeakIntensity~=0);
index2_Level1=find(msLEVEL==1);
index3_1=find(strcmp(ScanType,'Full')==1);
index4_1=intersect(index0,index2_Level1);
index_Level_1=intersect(index3_1,index4_1);
basePeaks=[basePeakMz(index_Level_1), basePeakIntensity(index_Level_1)];
Retention_Time_TIC1=Retention_Time(index_Level_1)/60;
Peaks_TIC1=cell(size(index_Level_1,1),1);
k=0;
for i=index_Level_1'
    k=k+1;
    Peaks_TIC1{k}=peaks{i};
end
save('Peaks_TIC1.mat','Peaks_TIC1','-v7.3')
save('Retention_Time_TIC1.mat','Retention_Time_TIC1')
save('basePeaks.mat','basePeaks')
%% Importing msLevel 2 peaks
index3_2=find(strcmp(ScanType,'Full')==1);
index2_Level_2=find(msLEVEL==2);
index4_2=intersect(index0,index2_Level_2);
index_Level_2=intersect(index3_2,index4_2);
Retention_Time_TIC2=Retention_Time(index_Level_2)/60;
Peaks_TIC2=cell(size(index_Level_2,1),1);
k=0;
for i=index_Level_2'
    k=k+1;
    Peaks_TIC2{k}=peaks{i};
end
save('Peaks_TIC2.mat','Peaks_TIC2','-v7.3')
save('Retention_Time_TIC2.mat','Retention_Time_TIC2')
%%
