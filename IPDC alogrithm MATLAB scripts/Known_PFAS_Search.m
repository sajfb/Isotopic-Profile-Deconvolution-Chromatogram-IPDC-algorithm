function Final_Candidates=Known_PFAS_Search(Final_Candidates,Column_xlx)
Known_PFCs_compounds={
    
% PerFluoroAlkyl Sulfonates (PFSAs)
'C16HF33O3S','TriTriaContaFluoroHexaDecaneSulfonic Acid (PFSA)'
'C15HF31O3S','HenTriaContaFluoroPentaDecaneSulfonic Acid (PFSA)'
'C14HF29O3S','NonaCosaFluoroTetraDecaneSulfonic Acid (PFSA)'
'C13HF27O3S','HeptaCosaFluoroTriDecaneSulfonic Acid (PFSA)'
'C12HF25O3S','PentaCosaFluoroDoDecaneSulfonic Acid (PFSA)'
'C11HF23O3S','TriCosaFluorounDecaneSulfonic Acid (PFSA)'
'C10HF21O3S','PerFluoroDecaneSulfonic Acid (PFSA)'
'C9HF19O3S','PerFluoroNonaneSulfonic Acid (PFSA)'
'C8HF17O3S','PerFluoroOctaneSulfonic Acid (PFSA)'
'C7HF15O3S','PerFluoroHeptaneSulfonic Acid (PFSA)'
'C6HF13O3S','PerFluoroHexaneSulfonic Acid (PFSA)'
'C5HF11O3S','PerFluoroPentaneSulfonic Acid (PFSA)'
'C4HF9O3S','PerFluoroButaneSulfonic Acid (PFSA)'

% PerFluoroAlkyl carboxylic Acids (PFCAs)
'C16HF31O2','PerFluoroPalmitic Acid (PFCA)'
'C15HF29O2','PerFluoroPentaDecanoic Acid (PFCA)'
'C14HF27O2','PerFluoroMyristic Acid (PFCA)'
'C13HF25O2','PerFluoroTriDecanoic Acid (PFCA)'
'C12HF23O2','PerFluoroDoDecanoic Acid (PFCA)'
'C11HF21O2','PerFluorouNDecanoic Acid (PFCA)'
'C10HF19O2','PerFluoroDecanoic Acid (PFCA)'
'C9HF17O2','PerFluoroNonanoic Acid (PFCA)'
'C8HF15O2','PerFluoroOctanoic Acid (PFCA)'
'C7HF13O2','PerFluoroHeptanoic Acid (PFCA)'
'C6HF11O2','PerFluoroHexanoic Acid (PFCA)'
'C5HF9O2','PerFluoroPentanoic Acid (PFCA)'
'C4HF7O2','PerFluoroButanoic Acid (PFCA)'

% PFCAs fragment [M-CHO2]-
'C15HF31','PFCAs fragment [M-CHO2]-'
'C14HF29','PFCAs fragment [M-CHO2]-'
'C13HF27','PFCAs fragment [M-CHO2]-'
'C12HF25','PFCAs fragment [M-CHO2]-'
'C11HF23','PFCAs fragment [M-CHO2]-'
'C10HF21','PFCAs fragment [M-CHO2]-'
'C9HF19','PFCAs fragment [M-CHO2]-'
'C8HF17','PFCAs fragment [M-CHO2]-'
'C7HF15','PFCAs fragment [M-CHO2]-'
'C6HF13','PFCAs fragment [M-CHO2]-'
'C5HF11','PFCAs fragment [M-CHO2]-'
'C4HF9','PFCAs fragment [M-CHO2]-'
'C3HF7','PFCAs fragment [M-CHO2]-'

% Mono-fluoroalkyl carboxylic acids (MFCAs)
'C16H31FO2','MonoFluoroPalmitic Acid (MFCA)'
'C15H29FO2','MonoFluoroPentaDecanoic Acid (MFCA)'
'C14H27FO2','MonoFluoroMyristic Acid (MFCA)'
'C13H25FO2','MonoFluoroTriDecanoic Acid (MFCA)'
'C12H23FO2','MonoFluoroDoDecanoic Acid (MFCA)'
'C11H21FO2','MonoFluoroNDecanoic Acid (MFCA)'
'C10H19FO2','MonoFluoroDecanoic Acid (MFCA)'
'C9H17FO2','MonoFluoroNonanoic Acid (MFCA)'
'C8H15FO2','MonoFluoroOctanoic Acid (MFCA)'
'C7H13FO2','MonoFluoroHeptanoic Acid (MFCA)'
'C6H11FO2','MonoFluoroHexanoic Acid (MFCA)'
'C5H9FO2','MonoFluoroPentanoic Acid (MFCA)'
'C4H7FO2','MonoFluoroButanoic Acid (MFCA)'

% Tri-fluoroalkyl carboxylic acids (TrFCAs)
'C16H29F3O2','TriFluoroPalmitic Acid (TrFCA)'
'C15H27F3O2','TriFluoroPentaDecanoic Acid (TrFCA)'
'C14H25F3O2','TriFluoroMyristic Acid (TrFCA)'
'C13H23F3O2','TriFluoroTriDecanoic Acid (TrFCA)'
'C12H21F3O2','TriFluoroDoDecanoic Acid (TrFCA)'
'C11H19F3O2','TriFluoroNDecanoic Acid (TrFCA)'
'C10H17F3O2','TriFluoroDecanoic Acid (TrFCA)'
'C9H15F3O2','TriFluoroNonanoic Acid (TrFCA)'
'C8H13F3O2','TriFluoroOctanoic Acid (TrFCA)'
'C7H11F3O2','TriFluoroHeptanoic Acid (TrFCA)'
'C6H9F3O2','TriFluoroHexanoic Acid (TrFCA)'
'C5H7F3O2','TriFluoroPentanoic Acid (TrFCA)'
'C4H5F3O2','TriFluoroButanoic Acid (TrFCA)'

% Tetra-fluoroalkyl carboxylic acids (TeFCAs)
'C10H16F4O2','TetraFluoroDecanoic Acid (TeFCA)'
'C9H14F4O2','TetraFluoroNonanoic Acid (TeFCA)'

% Penta-fluorodecanoic acid (PeFDA)
'C10H15F5O2','PentaFluoroDecanoic Acid (PeFDA)'

% Hexa- fluorodecanoic acid (HFCA)
'C10H14F6O2','HexaFluoroDecanoic Acid (HFCA)'

% Mono-fluoroalkyl ether carboxylic acids (MFECA)
'C5H9FO3','(MFECA)'
'C6H11FO3','(MFECA)'
'C7H13FO3','(MFECA)'
'C8H15FO3','(MFECA)'
'C9H17FO3','(MFECA)'
'C10H19FO3','(MFECA)'

% Tri-fluoroalkyl ether carboxylic acids (TrFECAs)
'C6H9F3O3','(TrFECAs)'
'C7H11F3O3','(TrFECAs)'
'C8H13F3O3','(TrFECAs)'
'C9H15F3O3','(TrFECAs)'
'C10H17F3O3','(TrFECAs)'

% Polyfluoroalkyl sulfonate
'C8H10F8O3S','OctaFluoroOctaneSulfonic Acid'
'C6H12F2O3S','DiFluoroHexaneSulfonic Acid'

% other known PFCs
'C8HF15O3S','Perfluoroethylcyclohexanesulfonate (PFECHS)'
'C4H2F9NO2S','perFluoro-1-Butane-SulfonAmide (FBSA)'
'C8H2F17NO2S','perFluoroOctaneSulfonAmide (FOSA)'
'C12H8F17NO4S','N-EthylperFluoroOctane SulfonAmide (N-EtFOSA)'
'C11H6F17NO4S','N-MethylperFluoroOctane SulfonamidoAcetic Acid (N-MeFOSAA)'
'C12H8F17NO4S','N-EthylperFluoroOctane SulfonamidoAcetic Acid (N-EtFOSAA)'
'C11H8F17NO3S','2(N-Methylperfluorooctanesulfonamido) Ethanol (N-MeFOSE)'
'C12H10F17NO3S','2(N-EthylperFluoroOctaneSulfonamido) Ethanol (N-EtFOSE)'

% PolyBromoPhenol (PBP)
'C6HBr5O','PentaBromoPhenol (PBP)'
'C6H2Br4O','TetraBromoPhenol (PBP)'
'C6H3Br3O','TriBromoPhenol (PBP)'
'C6H4Br2O','DiBromoPhenol (PBP)'
'C6H5BrO','MonoBromoPhenol (PBP)'

% PolyChloroPhenol (PCP)
'C6HCl5O','PentaChloroPhenol (PCP)'
'C6H2Cl4O','TetraChloroPhenol (PCP)'
'C6H3Cl3O','TriChloroPhenol (PCP)'
'C6H4Cl2O','DiChloroPhenol (PCP)'
'C6H5ClO','MonoChloroPhenol (PCP)'

% OH - PolyChlorinated Biphenyl (OH-PCB)
'C12Cl10O','OH - PerChlorinated Biphenyl (OH-PCB)'
'C12HCl9O','OH - NonaChlorinated Biphenyl (OH-PCB)'
'C12H2Cl8O','OH - OctaChlorinated Biphenyl (OH-PCB)'
'C12H3Cl7O','OH - HeptaChlorinated Biphenyl (OH-PCB)'
'C12H4Cl6O','OH - HexaChlorinated Biphenyl (OH-PCB)'
'C12H5Cl5O','OH - PentaChlorinated Biphenyl (OH-PCB)'
'C12H6Cl4O','OH - TetraChlorinated Biphenyl (OH-PCB)'
'C12H7Cl3O','OH - TriChlorinated Biphenyl (OH-PCB)'
'C12H8Cl2O','OH - DiChlorinated Biphenyl (OH-PCB)'
'C12H9ClO','OH - MonoChlorinated Biphenyl (OH-PCB)'

% OH - PolyBrominated Diphenyl Ethers (OH-BDE)
'C12Br10O2','OH - PerBrominated Diethyl Ether (OH-BDE)'
'C12HBr9O2','OH - NonaBrominated Diphenyl Ethers (OH-BDE)'
'C12H2Br8O2','OH - OctaBrominated Diphenyl Ethers (OH-BDE)'
'C12H3Br7O2','OH - HeptaBrominated Diphenyl Ethers (OH-BDE)'
'C12H4Br6O2','OH - HexaBrominated Diphenyl Ethers (OH-BDE)'
'C12H5Br5O2','OH - PentaBrominated Diphenyl Ethers (OH-BDE)'
'C12H6Br4O2','OH - TetraBrominated Diphenyl Ethers (OH-BDE)'
'C12H7Br3O2','OH - TriBrominated Diphenyl Ethers (OH-BDE)'
'C12H8Br2O2','OH - DiBrominated Diphenyl Ethers (OH-BDE)'
'C12H9BrO2','OH - MonoBrominated Diphenyl Ethers (OH-BDE)'

% OH - Polychlorinated Diethyl Ether (OH-CDE)
'C12Cl10O2','OH - PerChlorinated Diethyl Ether (OH-CDE)'
'C12HCl9O2','OH - NonaChlorinated Diethyl Ether (OH-CDE)'
'C12H2Cl8O2','OH - OctaChlorinated Diethyl Ether (OH-CDE)'
'C12H3Cl7O2','OH - HeptaChlorinated Diethyl Ether (OH-CDE)'
'C12H4Cl6O2','OH - HexaChlorinated Diethyl Ether (OH-CDE)'
'C12H5Cl5O2','OH - PentaChlorinated Diethyl Ether (OH-CDE)'
'C12H6Cl4O2','OH - TetraChlorinated Diethyl Ether (OH-CDE)'
'C12H7Cl3O2','OH - TriChlorinated Diethyl Ether (OH-CDE)'
'C12H8Cl2O2','OH - DiChlorinated Diethyl Ether (OH-CDE)'
'C12H9ClO2','OH - MonoChlorinated Diethyl Ether (OH-CDE)'

% OH- PolyChlorinated Diethyl Furan (OH-CDF)
'C12Cl8O2','OH- OctaChlorinated Diethyl Furan (OH-CDF)'
'C12HCl7O2','OH- HeptaChlorinated Diethyl Furan (OH-CDF)'
'C12H2Cl6O2','OH- HexaChlorinated Diethyl Furan (OH-CDF)'
'C12H3Cl5O2','OH- PentaChlorinated Diethyl Furan (OH-CDF)'
'C12H4Cl4O2','OH- TetraChlorinated Diethyl Furan (OH-CDF)'
'C12H5Cl3O2','OH- TriChlorinated Diethyl Furan (OH-CDF)'
'C12H6Cl2O2','OH- DiChlorinated Diethyl Furan (OH-CDF)'
'C12H7ClO2','OH- MonoChlorinated Diethyl Furan (OH-CDF)'

% Other Known Compounds
'C12H18Br6','Hexabromocyclododecane (HBCDD)'
'C15H12Br4O2','Tetrabromobisphenol A (TBBPA)'
'C7H3Br2NO','Bromoxynil'
'C9H6Cl6O3S','Endosulfan'
'C9H6Cl6O4S','Endosulfan Sulfate'
'C8HCl3N2O','OH-chlorothalonil'
'C8H2Cl2N2O','OH-chlorothalonil-b'
'C7H4Cl2O4','Dichloro dihydroxybenzoic acid (DCBA)'
'C7H4Br2O4','Dibromo dihydroxybenzoic acid (DBBA)'
};
%%
N_F=length(Column_xlx);
for n_f=1:N_F     % n_f represent the number of molecular adducts
    for i=1:size(Final_Candidates,1)
        i
        if ~isempty(Final_Candidates{i,Column_xlx(n_f)})
            x_Known_compounds=find(strcmp(Final_Candidates(i,Column_xlx(n_f)),Known_PFCs_compounds(:,1))==1);
            if ~isempty(x_Known_compounds)
                Final_Candidates{i,Column_xlx(n_f)+3}=Known_PFCs_compounds{x_Known_compounds,2};
            end
        end
    end
end
