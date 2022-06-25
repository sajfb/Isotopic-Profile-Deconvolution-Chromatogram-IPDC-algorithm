%% This module is an add-on to "ChromatographyOutput_..." to generate chromatogram trace for ID compounds
clc
ID=541076   %<-- ID of the compound of Interest
Mass_window=0.01;
Chromatogram_scan_interval=10;
Smoothing_window=10;
noise_removal=1;
Max_No_isomers=25;
%%
% load Final_Matches.mat
% load('Retention_Time_TIC1.mat')
% load('\Peaks_TIC1.mat')
% load ID_library.mat
%%
% INT=Final_Matches(:,4);
% ScanNumber=Final_Matches(:,7);
% EntireIDs=Final_Matches(:,1);
%%
x=find(ID==EntireIDs);
figure(ID)
MZ=Final_Matches(x(1),3); %Exact Mass
[PA, No_Isomers]=PeakAreaPlot(Peaks_TIC1,Retention_Time_TIC1,MZ,Mass_window,INT(x),ScanNumber(x),...
    Chromatogram_scan_interval,Smoothing_window,noise_removal,Max_No_isomers);
M=['m/z= ',num2str(MZ)];
annotation('textbox', [0.132651843419506 0.926743424113143 0.200352520818716 0.0643478246875432],...
    'String', M,'LineStyle','none','FontSize',16);
title({[num2str(ID),' - ','[',Chemical_structure_Print(ID_library(ID,:)),]})
