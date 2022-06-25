%% Molecular_Formula=[Carbon Hydrogen Bromine Chlorine Fluorine Nitrogen Sodium Oxygen Phosphorus Sulfur]
function MW_1=Most_Abundant_Mass(Molecular_Formula,X)
H=[1.00783 2.01410];
F=18.9984032;
I=126.90448;
N=[14.00307 15.00011];
Na=22.98976928;
O=[15.99491 16.99913 17.99916];
P=30.97377;
MW_1=Molecular_Formula(2)*H(1)+Molecular_Formula(5)*F...
    +Molecular_Formula(6)*I(1)+Molecular_Formula(7)*N(1)+...
    Molecular_Formula(8)*Na+Molecular_Formula(9)*O(1)+Molecular_Formula(10)*P+X;
