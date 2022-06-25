% x=find(5==ID_library(:,1)...    %C
%        & 4==ID_library(:,2)...   %H
%        & 0==ID_library(:,3)...   %Br
%        & 3==ID_library(:,4)...   %Cl
%        & 0==ID_library(:,5)...   %F
%        & 0==ID_library(:,6)...   %I
%        & 2==ID_library(:,7)...   %N
%        & 0==ID_library(:,8)...   %Na
%        & 1==ID_library(:,9)...   %O
%        & 0==ID_library(:,10)...  %P
%        & 1==ID_library(:,11))    %S
% Compound=ID_library(x,:);
function Chemical_structure=Chemical_structure_Print(Compound)
C=[];
if Compound(1)==1
    C='C';
elseif Compound(1)>1
    C=['C',num2str(Compound(1))];
end
H=[];
if Compound(2)==1
    H='H';
elseif Compound(2)>1
    H=['H',num2str(Compound(2))];
end
Br=[];
if Compound(3)==1
    Br='Br';
elseif Compound(3)>1
    Br=['Br',num2str(Compound(3))];
end
Cl=[];
if Compound(4)==1
    Cl='Cl';
elseif Compound(4)>1
    Cl=['Cl',num2str(Compound(4))];
end
F=[];
if Compound(5)==1
    F='F';
elseif Compound(5)>1
    F=['F',num2str(Compound(5))];
end
I=[];
if Compound(6)==1
    I='I';
elseif Compound(6)>1
    I=['I',num2str(Compound(6))];
end
N=[];
if Compound(7)==1
    N='N';
elseif Compound(7)>1
    N=['N',num2str(Compound(7))];
end
Na=[];
if Compound(8)==1
    Na='Na';
elseif Compound(8)>1
    Na=['Na',num2str(Compound(8))];
end
O=[];
if Compound(9)==1
    O='O';
elseif Compound(9)>1
    O=['O',num2str(Compound(9))];
end
P=[];
if Compound(10)==1
    P='P';
elseif Compound(10)>1
    P=['P',num2str(Compound(10))];
end
S=[];
if Compound(11)==1
    S='S';
elseif Compound(11)>1
    S=['S',num2str(Compound(11))];
end
Chemical_structure=[C,H,Br,Cl,F,I,N,Na,O,P,S]; % Hill notation system
% fprintf('%s \n',Chemical_structure)
