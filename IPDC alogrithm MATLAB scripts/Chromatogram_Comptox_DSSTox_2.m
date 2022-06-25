%% To create a folder containing of compounds that have at least one record in comptox for two ionization pathways
% clear;clc
% load Final_Matches.mat
% load Final_Candidates.mat
% load Retention_Time_TIC1.mat
% load('Peaks_TIC1.mat')
% load('XIC_Primary.mat')
ID_I=cell2mat(XIC_Primary(:,4));
ID_E=cell2mat(XIC_Primary(:,5));
%% these parameters should be exactly similar to Chromatographic analysis in ChromatographyOutput_...
Mass_window=0.01;
Chromatogram_scan_interval=10;
Smoothing_window=10;
%%
INT=Final_Matches(:,4);
ScanNumber=Final_Matches(:,7);
EntireIDs=Final_Matches(:,1);
mkdir('Chromatograms_Comptox_DSSTox')
currentFolder = pwd;
route=[currentFolder,'\Chromatograms_Comptox_DSSTox\'];
Column_xlx=[15, 19]; % Columns where the molecular formulas were placed
for i=1:size(Final_Candidates,1)
    if ~isempty(Final_Candidates{i,Column_xlx(1)}) || ~isempty(Final_Candidates{i,Column_xlx(2)})
        if ~isempty(Final_Candidates{i,3})     % The main isomer
            if Final_Candidates{i,11}>=2*Final_Candidates{i,3}+1     % More than 1 histogram
                EPA_ava=0;
                if Final_Candidates{i,Column_xlx(1)+2}>0
                    EPA_ava=1;
                elseif ~isempty(Final_Candidates{i,Column_xlx(2)+2})
                    if Final_Candidates{i,Column_xlx(2)+2}>0
                        EPA_ava=1;
                    end
                end
                Chromatogram_Shape=0;
                if Final_Candidates{i,3}<=5 % up to 5 isomers
                    Chromatogram_Shape=1;
                else % Or a known compound
                    if ~isempty(Final_Candidates{i,Column_xlx(1)+3}) || ~isempty(Final_Candidates{i,Column_xlx(2)+3})
                        Chromatogram_Shape=1;
                    end
                end
                if EPA_ava==1 && Chromatogram_Shape==1
                    ID=Final_Candidates{i,1}
                    x=find(EntireIDs==ID);
                    figure(ID)
                    axes1 = axes('FontSize',20);
                    x_N_Iso=find(ID>=ID_I & ID<=ID_E);
                    N_Iso=XIC_Primary{x_N_Iso,2};       % Number of isotopologues
                    if N_Iso>3
                        noise_removal=1;
                        Max_No_isomers=25;
                    else
                        noise_removal=10;
                        Max_No_isomers=5;
                    end
                    if ~isempty(Final_Candidates{i,Column_xlx(1)+3})
                        title({[num2str(i+1),' - ',num2str(ID),' - ',Final_Candidates{i,2}],Final_Candidates{i,Column_xlx(1)+3}})
                    elseif ~isempty(Final_Candidates{i,Column_xlx(2)+3})
                        title({[num2str(i+1),' - ',num2str(ID),' - ',Final_Candidates{i,2}],Final_Candidates{i,Column_xlx(2)+3}})
                    else
                        title({[num2str(i+1),' - ',num2str(ID),' - ',Final_Candidates{i,2}]});
                    end
                    MZ=Final_Matches(x(1),3); %Exact Mass
                    PeakAreaPlot(Peaks_TIC1,Retention_Time_TIC1,MZ,Mass_window,INT(x),ScanNumber(x),...
                        Chromatogram_scan_interval,Smoothing_window,noise_removal,Max_No_isomers);
                    M=['m/z= ',num2str(MZ)];
                    annotation('textbox', [0.132651843419506 0.926743424113143 0.200352520818716 0.0643478246875432],...
                        'String', M,'LineStyle','none','FontSize',16);
                    figname=[route,num2str(i+1),' - ',num2str(ID),' - ',Final_Candidates{i,2},'.png'];
                    pause(0.001);
                    frame_h = get(handle(gcf),'JavaFrame');
                    set(frame_h,'Maximized',1);
                    saveas(gcf,figname)
                    close(figure(ID))
                end
            end
        end
    end
end
