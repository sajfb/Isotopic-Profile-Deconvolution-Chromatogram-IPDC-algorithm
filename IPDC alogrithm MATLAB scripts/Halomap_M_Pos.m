load('Reduced_data.mat')
EntireKnownMolFormula={{'C12Cl10';'C12HCl9';'C12H2Cl8';'C12H3Cl7';'C12H4Cl6';'C12H5Cl5';'C12H6Cl4';'C12H7Cl3';'C12H8Cl2';'C12H9Cl'},'PCB',[1 0 0],'square';
    {'C12Cl10O';'C12HCl9O';'C12H2Cl8O';'C12H3Cl7O';'C12H4Cl6O';'C12H5Cl5O';'C12H6Cl4O';'C12H7Cl3O';'C12H8Cl2O';'C12H9ClO'},'PCDE',[1 0 1],'square';
    {'C12Cl8O';'C12HCl7O';'C12H2Cl6O';'C12H3Cl5O';'C12H4Cl4O';'C12H5Cl3O';'C12H6Cl2O';'C12H6ClO'},'PCDF',[0 0 1],'square';
    {'C12Br10O';'C12HBr9O';'C12H2Br8O';'C12H3Br7O';'C12H4Br6O';'C12H5Br5O';'C12H6Br4O';'C12H7Br3O';'C12H8Br2O';'C12H9BrO'},'PBDE','cyan','o';
    {'C8Cl8';'C8HCl7';'C8H2Cl6'},'PCStyrene',[0 0 1],'^';
    {'C6Cl6';'C6HCl5';'C6H2Cl4';'C6H3Cl3';'C6H4Cl2';'C6H5Cl'},'PCBz',[0.600000023841858 0.200000002980232 0],'v';
    {'C7H3Cl5O';'C7H4Cl4O';'C7H5Cl3O'},'PCA',[0.749019622802734 0.749019622802734 0],'o';
    {'C10Cl8';'C10HCl7';'C10H2Cl6';'C10H3Cl5';'C10H4Cl4';'C10H5Cl3';'C10H6Cl2';'C10H7Cl'},'PCN',[1 0 1],'o';
    {'C18H6Cl8';'C18H7Cl7';'C18H8Cl6';'C18H9Cl5';'C18H10Cl4';'C18H11Cl3';'C18H12Cl2';'C18H13Cl'},'PCT',[0 0 0],'x';
    {'C10H6Cl4O4'},'Dacthal',[0 0 0],'<';
    {'C14H9Cl5'},'DDT',[0 0 0],'o';
    {'C14H10Cl2O2'},'DDA',[0 0 0],'*';
    {'C14H8Cl4'},'DDE',[0 0.498039215803146 0],'^';
    {'C14H8Cl2'},'DDE [M-Cl_2]^+',[0 0.498039215803146 0],'diamond';
    {'C10H4Cl8O'},'Oxychlordane',[0 0 0],'pentagram';
    {'C10H5Cl7O'},'Oxychlordane [M+H-Cl]^+',[0 0 0],'pentagram';
    {'C10H5Cl8O'},'Oxychlordane [M+H]^+',[0 0 0],'pentagram';
    {'C14H9Cl3'},'DDMU',[0 0 0],'d';
    {'C10H5Cl7'},'Heptachlor',[0 0 0],'v';
    {'C10H6Cl8'},'Chlordane',[0 0 0],'>';
    {'C10H5Cl9'},'Nonachlor',[0 0 0],'h';
    {'C10H8Cl8'},'Toxaphene',[0 0 1],'*';
    {'C9H6Cl6O3S'},'Endosulfan',[1 0 0],'*';
    };
%% This section can be run once.
Column_xlx=[15 19]; % Columns where the molecular formulas were placed
KnownCompounds=cell(size(EntireKnownMolFormula,1),2);
K_H=[];
for k=1:size(EntireKnownMolFormula,1)
    k
    KnownMolFormula=EntireKnownMolFormula{k,1}
    for n_f=1:2     % n_f represent the number of ionization pathways
        for i=1:size(Final_Candidates,1)
            if ~isempty(Final_Candidates(i,Column_xlx(n_f)))
                x=find(strcmp(Final_Candidates(i,Column_xlx(n_f)),KnownMolFormula)==1, 1);
                if ~isempty(x)
                    if ~ischar(Final_Candidates{i,1})
                        ID_M=Final_Candidates{i,1};
                    end
                    No_Iso=Final_Candidates{i,3};
                    if No_Iso==1
                        if Final_Candidates{i,11}>2 && Final_Candidates{i,13}<0.8 && Final_Candidates{i,12}>25
                            K_H=[K_H;i,ID_M,Final_Candidates{i,4},Final_Candidates{i,6},Final_Candidates{i,9},Final_Candidates{i,10},Final_Candidates{i,14}];
                        end
                    else
                        for j=i+1:i+No_Iso
                            if Final_Candidates{j,11}>2 && Final_Candidates{j,13}<0.8 && Final_Candidates{j,12}>25
                                K_H=[K_H;j,ID_M,Final_Candidates{j,4},Final_Candidates{j,6},Final_Candidates{j,9},Final_Candidates{j,10},Final_Candidates{j,14}];
                            end
                        end
                    end
                end
            end
        end
        KnownCompounds{k,n_f}=K_H;
        K_H=[];
    end
end
save('KnownCompounds.mat','KnownCompounds')
%%
load KnownCompounds.mat
MplusOne=[0;1.9974;3.0007;3.9948];
Delta_mass=0.05;      % (Da)
Delta_RT=0.1;   % (min)
y=0;
G_common=[];G_common_i={};
for k=1:size(EntireKnownMolFormula,1)
    for n_f=1:2     % n_f represent the number of ionization pathways
        K_H=KnownCompounds{k,n_f};
        if ~isempty(K_H)
            i=1;
            while size(K_H,1)>=i
                x_MW=[];
                for M=1:length(MplusOne)
                    x_MW=[x_MW;find(abs(Reduced_data(:,3)-K_H(i,3)-MplusOne(M))<=Delta_mass)];
                    x_MW=[x_MW;find(abs(Reduced_data(:,3)-K_H(i,3)+MplusOne(M))<=Delta_mass)];
                end
                x_RT=find(abs(Reduced_data(:,5)-K_H(i,5))<=Delta_RT);
                x_MW_x_RT=intersect(x_MW,x_RT);
                if isempty(x_MW_x_RT)
                    K_H(i,:)=[];
                else
                    y=y+1;
                    G_common_i{y,1}=K_H(i,1);
                    G_common_i{y,2}=Reduced_data(x_MW_x_RT,1);
                    G_common=[G_common;Reduced_data(x_MW_x_RT,1)];
                    i=i+1;
                end
            end
            KnownCompounds{k,n_f}=K_H;
        end
    end
end
G_i=setdiff(Reduced_data(:,1),G_common);
L=[];
knownmasses=[];   % background masses due to PCB overloading
for i=1:size(Reduced_data,1)
    x=find(Reduced_data(i,1)==G_i, 1);
    if ~isempty(x)
        y=find(abs(Reduced_data(i,3)-knownmasses)<=Delta_mass, 1);
        if isempty(y)
            L=[L;Reduced_data(i,:)];
        end
    end
end
L2=[];
while size(L,1)~=0
    x_mz=find(abs(L(:,3)-L(1,3))<=Delta_mass);
    [~,x_P_area]=maxk(L(x_mz,4),5);    % Most abundant isomers for unknwon masses
    L2=[L2;L(x_mz(x_P_area),:)];
    L(x_mz,:)=[];
end
size(L2,1)
if~isempty(L2)
    L2=sortrows(L2,1);
else
    L2=[0 0 0 0 0 0 0];
end
figure(1)
axes1 = axes('FontSize',16);
box on
hold on
set(gca,'Color',[0.3010, 0.7450, 0.9330]);
plot(L2(:,5),L2(:,3),'DisplayName','Unknown features','MarkerFaceColor',[1 0.843137264251709 0],'MarkerEdgeColor',[1 0.843137264251709 0],'Marker','diamond',...
    'MarkerSize',4,'LineStyle','none',...
    'Color',[1 0.843137264251709 0]);
ylabel('m/z (Da)');
xlabel('Relative Retention Time (min)');
xlim([15 51.5]) % m/z limit
ylim([100 600]) % RT limit
set(gcf,'color','w');
Reduced_data_filtered=L2;
save('Reduced_data_filtered.mat','Reduced_data_filtered')
% title('')
%%
load('KnownCompounds.mat') % Optional command to see what known compounds did not pass the false positive prediction test
Number_known_features=0;
for k=1:size(EntireKnownMolFormula,1)
    for n_f=1:2     % n_f represent the number of ionization pathways
        K_H=KnownCompounds{k,n_f};
        Number_known_features=Number_known_features+size(K_H,1);
        if ~isempty(K_H)
            if n_f==1 && k~=14 && k~=17
                DP=[EntireKnownMolFormula{k,2},'  [M]^+'];
                plot(K_H(:,5),K_H(:,3),'DisplayName',DP,...
                    'MarkerFaceColor',EntireKnownMolFormula{k,3},'MarkerEdgeColor',EntireKnownMolFormula{k,3},'Marker',EntireKnownMolFormula{k,4},...
                    'MarkerSize',6,'LineStyle','none',...
                    'Color',EntireKnownMolFormula{k,3});
            elseif n_f==2 && k~=14 && k~=16
                DP=[EntireKnownMolFormula{k,2},'  [M-Cl]^+'];
                plot(K_H(:,5),K_H(:,3),'DisplayName',DP,...
                    'MarkerEdgeColor',EntireKnownMolFormula{k,3},'Marker',EntireKnownMolFormula{k,4},...
                    'MarkerSize',6,'LineStyle','none',...
                    'Color',EntireKnownMolFormula{k,3});
            else
                plot(K_H(:,5),K_H(:,3),'DisplayName',EntireKnownMolFormula{k,2},...
                    'MarkerFaceColor',EntireKnownMolFormula{k,3},'MarkerEdgeColor',EntireKnownMolFormula{k,3},'Marker',EntireKnownMolFormula{k,4},...
                    'MarkerSize',6,'LineStyle','none',...
                    'Color',EntireKnownMolFormula{k,3});
            end
        end
    end
end
legend1 = legend('show');
set(legend1,'Location','northwest','FontSize',12);
legend boxoff
%%
kn_str=[num2str(Number_known_features),' known features'];
annotation('textbox',...
    [0.710976027750288 0.140323564252958 0.242236018088293 0.0725126458963295],...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403],...
    'String',kn_str,...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');
%%
if L2(1,1)==0
else
    unkn_str=[num2str(size(L2,1)),' unknown features'];
    annotation('textbox',...
        [0.710976027750289 0.172261890244148 0.242236018088293 0.0725126458963295],...
        'Color',[0 0.498039215803146 0],...
        'String',unkn_str,...
        'LineStyle','none',...
        'FontSize',16,...
        'FitBoxToText','off');
end
%%
pause(0.000001)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
%%
datacursorextra % a tool to use datatips more efficiently on MATLAB figures
%%
h = findobj(gca,'Type','line');
y=get(h,'Ydata');y=y(end:-1:1);
known=0;
p=length(y);
for i=2:p
    known=known+length(find(~isnan(y{i})));
    y{i,2}=length(find(~isnan(y{i})));
end
known
unknown=length(find(~isnan(y{1})))
