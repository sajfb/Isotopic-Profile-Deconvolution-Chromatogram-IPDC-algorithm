load('Reduced_data.mat')
EntireKnownMolFormula={{'C12Cl10';'C12HCl9';'C12H2Cl8';'C12H3Cl7';'C12H4Cl6';'C12H5Cl5';'C12H6Cl4';'C12H7Cl3';'C12H8Cl2';'C12H9Cl'},'PCBs',[1 0 0],'square';
    {'C16HF33O3S';'C15HF31O3S';'C14HF29O3S';'C13HF27O3S';'C12HF25O3S';'C11HF23O3S';'C10HF21O3S';'C9HF19O3S';'C8HF17O3S';'C7HF15O3S';'C6HF13O3S';'C5HF11O3S';'C4HF9O3S'},'PFAS',[1 0 1],'square';
    {'C16HF31O2';'C15HF29O2';'C14HF27O2';'C13HF25O2';'C12HF23O2';'C11HF21O2';'C10HF19O2';'C9HF17O2';'C8HF15O2';'C7HF13O2';'C6HF11O2';'C5HF9O2';'C4HF7O2'},'PFAAs',[0 0 1],'square';
    {'C13HF27';'C12HF25';'C11HF23';'C10HF21';'C9HF19';'C8HF17';'C7HF15';'C6HF13';'C5HF11';'C4HF9'},'PFAAs [M-CO_2]^-','g','o';
    {'C16H31FO2';'C15H29FO2';'C14H27FO2';'C13H25FO2';'C12H23FO2';'C11H21FO2';'C10H19FO2';'C9H17FO2';'C8H15FO2';'C7H13FO2';'C6H11FO2';'C5H9FO2';'C4H7FO2'},'MFCA',[0 1 0],'square';
    {'C16H29F3O2';'C15H27F3O2';'C14H25F3O2';'C13H23F3O2';'C12H21F3O2';'C11H19F3O2';'C10H17F3O2';'C9H15F3O2';'C8H13F3O2';'C7H11F3O2';'C6H9F3O2';'C5H7F3O2';'C4H5F3O2'},'TrFCA',[0 0 0],'o';
    {'C10H16F4O2';'C9H14F4O2'},'TeFCA',[0 0.498039215803146 0],'o';
    {'C10H15F5O2'},'PeFDA','w','o';
    {'C10H14F6O2'},'HFCA',[0 0 1],'^';
    {'C5H9FO3';'C6H11FO3';'C7H13FO3';'C8H15FO3';'C9H17FO3';'C10H19FO3'},'MFECA',[0.600000023841858 0.200000002980232 0],'v';
    {'C6H9F3O3';'C7H11F3O3';'C8H13F3O3';'C9H15F3O3';'C10H17F3O3'},'TrFECAs',[0 0.498039215803146 0],'diamond';
    {'C6H12F2O3S'},'DiFluoroHexaneSulfonic Acid','w','>';
    {'C8H10F8O3S',},'OctaFluoroOctaneSulfonic Acid','w','h';
    {'C8HF15O3S',},'PFECHS',[0 0 1],'*';
    };
%% This section can be run once.
Column_xlx=15; % Columns where the molecular formulas were placed
KnownCompounds=cell(size(EntireKnownMolFormula,1),1);
load('C:\Users\ACD\Desktop\APGC\Noise_Neg_mz_RT_Blood.mat') % detected features in blank
Delta_mass=0.01;  % (Da)
Delta_RT=1; % (min)
K_H=[];
for k=1:size(EntireKnownMolFormula,1)
    KnownMolFormula=EntireKnownMolFormula{k,1}
    for n_f=1:1     % n_f represent the number of molecular adducts
        for i=1:size(Final_Candidates,1)
            if ~isempty(Final_Candidates(i,Column_xlx(n_f)))
                x=find(strcmp(Final_Candidates(i,Column_xlx(n_f)),KnownMolFormula)==1, 1);
                if ~isempty(x)
                    if ~ischar(Final_Candidates{i,1})
                        ID_M=Final_Candidates{i,1};
                    end
                    No_Iso=Final_Candidates{i,3};
                    if No_Iso==1
                        y=find(abs(Noise_Neg_mz_RT_Blood(:,1)-Final_Candidates{i,4})<=Delta_mass & abs(Noise_Neg_mz_RT_Blood(:,2)-Final_Candidates{i,9})<=Delta_RT, 1);
                        if isempty(y)
                            K_H=[K_H;i,ID_M,Final_Candidates{i,4},Final_Candidates{i,6},Final_Candidates{i,9},Final_Candidates{i,10},Final_Candidates{i,14}];
                        end
                    else
                        for j=i+1:i+No_Iso
                            y=find(abs(Noise_Neg_mz_RT_Blood(:,1)-Final_Candidates{j,4})<=Delta_mass & abs(Noise_Neg_mz_RT_Blood(:,2)-Final_Candidates{j,9})<=Delta_RT, 1);
                            if isempty(y)
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
load('KnownCompounds.mat')
IAcc=0.01;
delta_t=0.1;
y=0;
G_common=[];G_common_i={};
for k=1:size(EntireKnownMolFormula,1)
    K_H=KnownCompounds{k,1};
    if ~isempty(K_H)
        i=1;
        while size(K_H,1)>=i
            x_MW=find(abs(Reduced_data(:,3)-K_H(i,3))<=IAcc);
            x_RT=find(abs(Reduced_data(:,5)-K_H(i,5))<=delta_t);
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
        KnownCompounds{k,1}=K_H;
    end
end
G_i=setdiff(Reduced_data(:,1),G_common);
L=[];
knownmasses=[];
for i=1:size(Reduced_data,1)
    x=find(Reduced_data(i,1)==G_i);
    if ~isempty(x)
        y=find(abs(Reduced_data(i,3)-knownmasses)<=IAcc);
        if isempty(y)
            L=[L;Reduced_data(i,:)];
        end
    end
end
L2=[];
while size(L,1)~=0
    x_mz=find(abs(L(:,3)-L(1,3))<=IAcc);
    [~,x_P_area]=maxk(L(x_mz,4),7);    % Most abundant isomers
    L2=[L2;L(x_mz(x_P_area),:)];
    L(x_mz,:)=[];
end
size(L2,1)
figure(2)
if~isempty(L2)
    L2=sortrows(L2,1);
else
    L2=[0 0 0 0 0 0 0 0 0 0 0 0 0 ]
end
axes1 = axes('FontSize',16);
box on
plot(L2(:,5),L2(:,3),'DisplayName','Unknown features','MarkerFaceColor',[1 0.843137264251709 0],'MarkerEdgeColor',[1 0.843137264251709 0],'Marker','diamond',...
    'MarkerSize',4,'LineStyle','none',...
    'Color',[1 0.843137264251709 0]);
ylabel('m/z (Da)');
xlabel('RT (min)');
set(axes1,'Color',[0.301 0.745 0.933],'FontSize',16);
title('1ng of PFAA and PFAS standard')
%%
hold on
load('KnownCompounds.mat') % Optional command to see what known compounds did not pass the false positive prediction test
Number_known_features=0;
for k=1:size(EntireKnownMolFormula,1)
    K_H=KnownCompounds{k,1};
    Number_known_features=Number_known_features+size(K_H,1);
    if ~isempty(K_H)
        if k~=4
            DP=[EntireKnownMolFormula{k,2},'  [M-H]^-'];
        else
            DP=[EntireKnownMolFormula{k,2}];
        end
        plot(K_H(:,5),K_H(:,3),'DisplayName',DP,...
            'MarkerFaceColor',EntireKnownMolFormula{k,3},'MarkerEdgeColor',EntireKnownMolFormula{k,3},'Marker',EntireKnownMolFormula{k,4},...
            'MarkerSize',6,'LineStyle','none',...
            'Color',EntireKnownMolFormula{k,3});
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
datacursorextra
