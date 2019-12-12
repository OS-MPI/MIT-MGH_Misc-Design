function [CoilGeom] = ToroidOptimizer(Dia,L)
%Summary: To use this function enter the wire diameter and the target
%inductance (in meters and Henries). The calculations assume two fully
%wound layers as it is a good balance between heat dissipation and
%geometric efficiency. The output is a struct containing the geometric
%information for a toroid with a D shaped core, and a circular core. 


%%Constants
mu0 = 4*pi*10^-7;
L0 = mu0*Dia/(2*pi);
L = L/4; %Use the assumption that there will be two fully wound layers...
%If there are two fully wound layers, it can be calculated at just one
%layer with 1/4 the inductnace. Then double it when winding
LperL0 = L/L0;
CoilGeom.General.SingleLayerInductance = L/4;
CoilGeom.General.WireDiameter = Dia;
CoilGeom.General.TargetInductance = L;

%% D shaped core
Alpha = 3; %User can change this to an integer between 2 and 5 (inclusive) varies the aspect ratio of the D
%A high alpha value such as 5 will be a taller/bigger core with fewer
%turns. It is theoretically more efficient but will probably have more
%leakage flux. 

% AlphaMat Key. The columns represent [alpha e h p s t]

AlphaMat = [2 0.26 0.645 3.6  0.72 1.77
    3 0.85 1.5   8.0  2.74 4.42
    4 1.6  2.4   12.8 5.76 8.09
    5 2.45 3.4   17.9 9.61 12.6]; %Table 1 From  "D-Shaped toroidal cage inductors" P.N. Murgatroyd & D. Belahrache,1989

T = AlphaMat(Alpha-1,6); %Picking a value in the table for clarity
S = AlphaMat(Alpha-1,5); %Also Selecting value for clarity
P = AlphaMat(Alpha-1,4); % As above

KFunk_D = @(k) (2*pi)^.5*(S/(P^(3/2)))*k.^(3/2)+0.25*k-LperL0; %From "The Optimal Form for Coreless Inductor", P. Murgatroyd IEEE TMI, 25 No. 3 1989
[K,~] =mybisect(KFunk_D, 0,1e5,50); %solving for ideal K with bisection method
WireLength_D = K*Dia;
Turns_D = (2*pi*K/P)^(1/2); % Also from P. Murgatroyd, section 2 in "Economic designs for single-layer toroidal inductors"
% TurnsD = sqrt(2*pi*L/(T*mu0*B)); %Turns needed in a D shaped core toroid
% ID_DToroid = 2*B;
B = WireLength_D/(Turns_D*P);

%% For the following reference Fig 4 in "D-Shaped toroidal cage inductors" P.N. Murgatroyd & D. Belahrache,1989
CoilGeom.DCore.FlatHeight = B*AlphaMat(Alpha-1,2);
CoilGeom.DCore.MaxHeight = B*AlphaMat(Alpha-1,2);
CoilGeom.DCore.RadiusAtPeak = B+Alpha^.5*B;
CoilGeom.DCore.ID = B*2;
CoilGeom.DCore.OD = 2*Alpha*B;
CoilGeom.DCore.Turns = Turns_D*2;
CoilGeom.DCore.Layers = 2;
CoilGeom.DCore.WireLength = WireLength_D;


%% Two-layer circular core


KFunk = @(k) 0.2722*k.^(3/2)+0.25*k-LperL0; %From "The Optimal Form for Coreless Inductor", P. Murgatroyd IEEE TMI, 25 No. 3 1989
[K,~] =mybisect(KFunk, 0,1e5,50); %solving for ideal K with bisection method
WireLength = K*Dia;
Turns = 0.8165*K^(1/2); % Also from P. Murgatroyd, section 2 in "Economic designs for single-layer toroidal inductors"
R = WireLength/(2*pi*Turns);% Equation 2 in "Economic designs for single-layer toroidal inductors",P. Murgatroyd
T = Dia/(sin(pi/Turns)*2)+R;% Equation 3 in "Economic designs for single-layer toroidal inductors"

CoilGeom.Circle.Turns = Turns*2;
CoilGeom.Circle.WireLength = WireLength;
CoilGeom.Circle.Layers=2;
CoilGeom.Circle.ID = 2*(T-R);
CoilGeom.Circle.OD = 2*(T+R);
CoilGeom.Circle.CoreRadius = R;
CoilGeom.Circle.CenterRadius = T;


%% Back-checking the work with two approximations. Laa is from hyperphysics, Lbb is from the paper that these geometries were derived from

Laa = mu0*Turns^2*R^2/(2*T);%Another formula for inductance of circular
% toroid from hyperphysics
Lbb = mu0*Dia/(2*pi)*(Turns.^2*(Turns+K/Turns-sqrt(Turns^2+2*K))+K/4);

end

