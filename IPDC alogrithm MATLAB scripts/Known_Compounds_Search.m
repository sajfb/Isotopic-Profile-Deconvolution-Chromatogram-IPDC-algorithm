function Final_Candidates=Known_Compounds_Search(Final_Candidates,Column_xlx)
Known_compounds={
    
% PCBs
'C12Cl10','PerChloroBiphenyl (PCB)'
'C12HCl9','NonaChloroBiphenyl (PCB)'
'C12H2Cl8','OctaChloroBiphenyl (PCB)'
'C12H3Cl7','HeptaChloroBiphenyl (PCB)'
'C12H4Cl6','HexaChloroBiphenyl (PCB)'
'C12H5Cl5','PentaChloroBiphenyl (PCB)'
'C12H6Cl4','TetraChloroBiphenyl (PCB)'
'C12H7Cl3','TriChloroBiphenyl (PCB)'
'C12H8Cl2','DiChloroBiphenyl (PCB)'
'C12H9Cl','MonoChloroBiphenyl (PCB)'
'C12H10','Biphenyl (CB)'

% PCDEs
'C12Cl10O','PerChloroDiphenylEther (PCDE)'
'C12HCl9O','NonaChloroDiphenylEther (PCDE)'
'C12H2Cl8O','OctaChloroDiphenylEther (PCDE)'
'C12H3Cl7O','HeptaChloroDiphenylEther (PCDE)'
'C12H4Cl6O','HexaChloroDiphenylEther (PCDE)'
'C12H5Cl5O','PentaChloroDiphenylEther (PCDE)'
'C12H6Cl4O','TetraChloroDiphenylEther (PCDE)'
'C12H7Cl3O','TriChloroDiphenylEther (PCDE)'
'C12H8Cl2O','DiChloroDiphenylEther (PCDE)'
'C12H9ClO','MonoChloroDiphenylEther (PCDE)'

% PCDFs
'C12Cl8O','PerChloroDibenzoFuran (PCDF)'
'C12HCl7O','HeptaChloroDibenzoFuran (PCDF)'
'C12H2Cl6O','HexaChloroDibenzoFuran (PCDF)'
'C12H3Cl5O','PentaChloroDibenzoFuran (PCDF)'
'C12H4Cl4O','TetraChloroDibenzoFuran (PCDF)'
'C12H5Cl3O','TriChloroDibenzoFuran (PCDF)'
'C12H6Cl2O','DiChloroDibenzoFuran (PCDF)'
'C12H6ClO','MonoChloroDibenzoFuran (PCDF)'

% PCDDs
'C12Cl8O2','[M-H-Cl+O2]- fragments for PCBs'
'C12HCl7O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H2Cl6O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H3Cl5O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H4Cl4O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H5Cl3O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H6Cl2O2','[M-H-Cl+O2]- fragments for PCBs'
'C12H7ClO2','[M-H-Cl+O2]- fragments for PCBs'

%PBDEs
'C12Br10O','PerBromoDiphenylEther (PBDE)'
'C12HBr9O','NonaBromoDiphenylEther (PBDE)'
'C12H2Br8O','OctaBromoDiphenylEther (PBDE)'
'C12H3Br7O','HeptaBromoDiphenylEther (PBDE)'
'C12H4Br6O','HexaBromoDiphenylEther (PBDE)'
'C12H5Br5O','PentaBromoDiphenylEther (PBDE)'
'C12H6Br4O','TetraBromoDiphenylEther (PBDE)'
'C12H7Br3O','TriBromoDiphenylEther (PBDE)'
'C12H8Br2O','DiBromoDiphenylEther (PBDE)'
'C12H9BrO','MonoBromoDiphenylEther (PBDE)'

% PBDFs
'C12Br8O','[M-H-Br2+O]- fragments for PBDEs'
'C12HBr7O','[M-H-Br2+O]- fragments for PBDEs'
'C12H2Br6O','[M-H-Br2+O]- fragments for PBDEs'
'C12H3Br5O','[M-H-Br2+O]- fragments for PBDEs'
'C12H4Br4O','[M-H-Br2+O]- fragments for PBDEs'
'C12H5Br3O','[M-H-Br2+O]- fragments for PBDEs'
'C12H6Br2O','[M-H-Br2+O]- fragments for PBDEs'
'C12H7BrO','[M-H-Br2+O]- fragments for PBDEs'

% Poly and per-chloro styrenes
'C8Cl8','PerChloroStyrene'
'C8HCl7','HeptaChloroStyrene'
'C8H2Cl6','HexaChloroStyrene'
'C8H6Cl2','DiChloroStyrene'
'C10H7Cl7','Dihydroheptachlor'

% Halogenated benzenes
'C6Cl6','HexaChloroBenzene'
'C6HCl5','PentaChloroBenzene'
'C6H2Cl4','TetraChloroBenzene'
'C6H3Cl3','TriChloroBenzene'
'C6H4Cl2','DiChloroBenzene'
'C6H5Cl','MonoChloroBenzene'
'C6H3Br3','TriBromoBenzene'

% PolyChloro anthracene
'C14H7Cl3','TriChloroAnthracene'
'C14H8Cl2','DiChloroAnthracene'

% PolyHalo anisole
'C7H3Cl5O','PentaChloroAnisole'
'C7H4Cl4O','TetraChloroAnisole'
'C7H5Cl3O','TriChloroAnisole'
'C7H5Br3O','TriBromoAnisole'
'C7H6Br2O','DiBromoAnisole'

% PolyChloro naphthalene
'C10Cl8','PerChloroNaphthalene'
'C10HCl7','HeptaChloroNaphthalene'
'C10H2Cl6','HexaChloroNaphthalene'
'C10H3Cl5','PentaChloroNaphthalene'
'C10H4Cl4','TetraChloroNaphthalene'
'C10H5Cl3','TriChloroNaphthalene'
'C10H6Cl2','DiChloroNaphthalene'
'C10H7Cl','MonoChloroNaphthalene'

% PolyChloroTerphenyl
'C18H5Cl9','NonaChloroTerphenyl'
'C18H6Cl8','OctaChloroTerphenyl'
'C18H7Cl7','HeptaChloroTerphenyl'
'C18H8Cl6','HexaChloroTerphenyl'
'C18H9Cl5','PentaChloroTerphenyl'
'C18H10Cl4','TetraChloroTerphenyl'
'C18H11Cl3','TriChloroTerphenyl'
'C18H12Cl2','DiChloroTerphenyl'
'C18H13Cl','MonoChloroTerphenyl'

% PolyChloroThioAnisole
'C7H3Cl5S','PentaChloroThioAnisole'
'C7H4Cl4S','TetraChloroThioAnisole'
'C7H5Cl3S','TriChloroThioAnisole'
'C7H6Cl2S','DiChloroThioAnisole'
'C7H7ClS','MonoChloroThioAnisole'

% Chrlodane derivatives
'C10H5Cl7','Heptachlor'
'C10H6Cl8','Chlordane'
'C10H5Cl9','Nonachlor'
'C10H4Cl8O','Oxychlordane'
'C12H9Cl7','Heptachlorotetracyclododecene'
'C10H8Cl8','Toxaphene'

% Other Known Compounds
'C10H6Cl4O4','Dacthal'
'C14H9Cl5','DDT'
'C14H10Cl2O2','DDA'
'C14H8Cl4','DDE'
'C14H10Cl4','DDD'
'C10H5Cl7O','Heptachlor Epoxide'
'C14H9Cl3','DDMU'
'C12H8Cl6O','Dieldrin'
'C18H10Cl2','DiChloroTetracene'
'C7H5Br3','TriBromoToluene'
'C9H6Cl6O3S','Endosulfan'
'C9H6Cl6O4S','Endosulfan Sulfate'
'C6H3Br3O','Tribromophenol'
'C14H9Cl5O','Dicofol'
'C6HCl5O','Pentachlorophenol'
'C9H6Cl6O3S','Endosulfan'


% Some identified Compounds
'C7H6Cl2O2','Di-chloromethoxyphenol'
'C7H7ClO2','Mono-chloromethoxyphenol'
'C7H7BrO2','Mono-bromomethoxyphenol'
'C9H9ClO3','Mono-chloromethoxyphenol acetate'
'C11H15ClO','Mono-chlordimethylpropylphenol'
'C9H11Cl2O3P','Probably a Metabolite'
'C12H9Cl2NO3','Vinclozoline'
'C12H7Cl2NO3','Nitrofen'
'C12H8Br2SO2','Bromophenyl Sulfone'
'C12H8Cl2SO2','Bis(4-chlorophenyl) sulfone'
'C8H16Cl3O4P','Probably a Metabolite'
'C9H10Cl2N2O','Diuron'
'C13H4Cl8','A Common Mass'
'C13H6Cl6','A Common Mass'
'C13H11Br2N','A Common Mass'
'C13H8Cl2O','DichloroBenzophenone'
'C18H10Cl4','TetraChloroTerphenyl'
'C10H6Cl2N2','Fenclorim'
'C14H18Cl2O4','Esteron 99'
'C10Cl10','Dienochlor'
'C9H2Cl4O2','TCID'
'C7HCl4N','TetraChloroBenzoNitrile'
'C13H6Cl2O2','Dichloroxanthen-9-one'
'C7H3Cl5','Dichlorobenzotrichloride'
'C7H7Cl2SNO2','Dichloramine T'
'C13H6Cl2O2','Dichloroxanthen-9-one'
'C23H18Cl2N2O2','TCMDC-123477'
'C28H16Cl2','Dichloro-bianthracene'
'C12H11Cl3O3','2,4-D chlorocrotyl ester'
'C8Cl6','Pentachloro-6-(chloroethynyl)benzene'
'C7H3Cl3O2','Trichlorobenzoic acid'
'C18H13Cl2N3','Climazolam'
'C9H11Cl3NO3PS','Chlorpyrifos'
};
%%
for n_f=1:length(Column_xlx)     % n_f represent the number of molecular adducts
    for i=1:size(Final_Candidates,1)
        i
        if ~isempty(Final_Candidates{i,Column_xlx(n_f)})
            x_Known_compounds=find(strcmp(Final_Candidates(i,Column_xlx(n_f)),Known_compounds(:,1))==1);
            if ~isempty(x_Known_compounds)
                Final_Candidates{i,Column_xlx(n_f)+3}=Known_compounds{x_Known_compounds,2};
            end
        end
    end
end
