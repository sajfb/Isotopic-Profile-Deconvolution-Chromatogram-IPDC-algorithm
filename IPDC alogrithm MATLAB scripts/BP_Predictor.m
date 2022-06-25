%% The function to predict the compound boiling point
function Tb=BP_Predictor(m,cl,br)
A=[52.54534071466831	0.47117954907202964	17.15965336330145	1.0718218519565328	80.5543947221509	0.9491977095603943];
a=A(1);b=A(2);c=A(3);d=A(4);e=A(5);f=A(6);
Tb=a*(m-cl*35.453-br*79.904).^b+c*cl.^d+e*br.^f-273.15;
end
