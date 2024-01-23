% ENME302
% Thomas Mitchell
% General case bar

clear
clc
close all

% ------------------------------ Parameters -------------------------------
Disp_mag = 1000; % Magnification factor
N_points = 10; % Number of points on plot
E = 200*10^9; % Young's modulus in Pa
A = pi*0.1^2/4; % Cross-sectional area of bar in m^2
Q_nodal = [0; 100000]; % Global forcing terms in N

% Element lengths in m:
L1 = 10;
L2 = 1.41*L1;

% Transformations angles in degrees:
alpha1 = 0;
alpha2 = 45;

% Assembly matrices:
AM1 = [0 0 1 0
       0 0 0 1];

AM2 = [0 0 1 0
       0 0 0 1];

% ---------------------------- Stiffness Terms ----------------------------
% Element 1:
K1 = local_bar(E,A,L1); % Element stiffness matrix in local coords
[K1hat, Lambda1] = global_bar(K1,alpha1); % Global coords
KG1 = assembly_bar(K1hat,AM1); % Structural contribution

% Element 2:
K2 = local_bar(E,A,L2); % Element stiffness matrix in local coords
[K2hat, Lambda2] = global_bar(K2,alpha2); % Global coords
KG2 = assembly_bar(K2hat,AM2); % Structural contribution

% --------------------------- Structural Terms ----------------------------
KG = KG1+KG2; % Structural stiffness matrix
q = KG\Q_nodal; % Structural displacements 

% --------------------- Displacements and Reactions -----------------------
% Element 1:
D1 = AM1'*q; % Global deflections
d1 = Lambda1*D1; % Local deflections
f1 = K1*d1; % Local nodal forces
F1 = K1hat*D1; % Global nodal forces

% Element 2:
D2 = AM2'*q; % Global deflections
d2 = Lambda2*D2; % Local deflections
f2 = K2*d2; % Local nodal forces
F2 = K2hat*D2; % Global nodal forces

% --------------------------- Stress and Strain ---------------------------
strain1 = (d1(2)-d1(1))/L1;
strain2 = (d2(2)-d2(1))/L2;

stress1 = f1(1)/A;
stress2 = f2(1)/A;

% --------------------------------- Plot ----------------------------------
% original structure
line([0 10],[0 10])
line([0 10],[10 10])

% structure with deformation (magnified)
line([0 10+Disp_mag*q(1)],[0 10+Disp_mag*q(2)])
line([0 10+Disp_mag*q(1)],[10 10+Disp_mag*q(2)])