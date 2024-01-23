% ENME302
% Thomas Mitchell
% General case frame

clear
clc
close all

% ------------------------------ Parameters -------------------------------
Disp_mag = 20; % Magnification factor
N_points = 10; % Number of points on plot
E = 200*10^9; % Young's modulus in Pa
I = 7.5*10^-6; % 2nd moment of area of frame
A = 9*10^-4; % Cross-sectional area of frame in m^2
Q_nodal = [0; 0; 0; 0; 0; 0]; % Global forcing terms in N

% Element lengths in m:
L1 = 2.5;
L2 = 2.5;
L3 = 1;

% Transformations angles in degrees:
alpha1 = -90;
alpha2 = -36.87;
alpha3 = 90;

% Assembly matrices:
AM1 = [1 0 0 0 0 0; 
       0 1 0 0 0 0; 
       0 0 1 0 0 0; 
       0 0 0 0 0 0; 
       0 0 0 0 0 0; 
       0 0 0 0 0 0];

AM2 = [1 0 0 0 0 0; 
       0 1 0 0 0 0; 
       0 0 1 0 0 0; 
       0 0 0 1 0 0; 
       0 0 0 0 1 0; 
       0 0 0 0 0 1];

AM3 = [0 0 0 0 0 0; 
       0 0 0 0 0 0; 
       0 0 0 0 0 0; 
       0 0 0 1 0 0; 
       0 0 0 0 1 0; 
       0 0 0 0 0 1];

% Complex loads:
wbar_UDL = 0; % UDL intensity
wbar_LVL = 1000*9.81*3*2.5; % LVL peak intensity
wbar_TPL = 0; % Transverse point load
pbar_DAL = 0; % Distributed axial load
pbar_CAL = 0; % Concentrated axial load

% Equivalent nodal vectors
Q_UDL = zeros(6,1);
Q_LVL = zeros(6,1);
Q_TPL = zeros(6,1);
Q_DAL = zeros(6,1);
Q_CAL = zeros(6,1);

% ---------------------------- Stiffness Terms ----------------------------
% Element 1:
K1 = local_frame(E,I,A,L1); % Element stiffness matrix in local coords
[K1hat, Lambda1] = global_frame(K1,alpha1); % Global coords
KG1 = assembly_frame(K1hat,AM1); % Structural contribution

% Element 2:
K2 = local_frame(E,I,A,L2); % Element stiffness matrix in local coords
[K2hat, Lambda2] = global_frame(K2,alpha2); % Global coords
KG2 = assembly_frame(K2hat,AM2); % Structural contribution

% Element 3:
K3 = local_frame(E,I,A,L3); % Element stiffness matrix in local coords
[K3hat, Lambda3] = global_frame(K3,alpha3); % Global coords
KG3 = assembly_frame(K3hat,AM3); % Structural contribution

% --------------------------- Equivalent Loads ----------------------------
% Uniformly distributed load:
% [f_UDL,F_UDL,Q_UDL] = equiv_uniform_load(wbar_UDL,L1,Lambda1,AM1);

% Linearly varying load:
[f_LVL,F_LVL,Q_LVL] = equiv_linear_load(wbar_LVL,L1,Lambda1,AM1);

% Transverse point load:
% [f_TPL,F_TPL,Q_TPL] = equiv_point_load(wbar_TPL,L1,L1/2,Lambda1,AM1);

% Distributed axial load:
% [f_DAL,F_DAL,Q_DAL] = equiv_distributed_axial_load(pbar_DAL,L1,Lambda1,AM1);

% Concentrated axial load:
% [f_CAL,F_CAL,Q_CAL] = equiv_concentrated_axial_load(pbar_CAL,L1,L1/2,Lambda1,AM1);

% --------------------------- Structural Terms ----------------------------
KG = KG1+KG2+KG3; % Structural stiffness matrix
Q_total = Q_nodal+Q_UDL+Q_LVL+Q_TPL+Q_DAL+Q_CAL; % Structural forces
q = KG\Q_total; % Structural displacements

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

% Element 3:
D3 = AM3'*q; % Global deflections
d3 = Lambda3*D3; % Local deflections
f3 = K3*d3; % Local nodal forces
F3 = K3hat*D3; % Global nodal forces

% --------------------------- Stress and Strain ---------------------------
strain1 = (d1(4)-d1(1))/L1;
strain2 = (d2(4)-d2(1))/L2;
strain3 = (d3(4)-d3(1))/L3;

stress1 = f1(1)/A;
stress2 = f2(1)/A;
stress3 = f3(1)/A;

% Axial and transverse displacements in local and global coords
[d_x,d_y,D_x,D_y] = deflections(d2,L2,L2/2,alpha2);

% --------------------------------- Plot ----------------------------------
hold on
Plot_deflected_shape(0,L1,0,0,d1,L1,alpha1,Disp_mag,N_points);
Plot_deflected_shape(0,L1,2,L3,d2,L2,alpha2,Disp_mag,N_points);
Plot_deflected_shape(2,0,2,L3,d3,L3,alpha3,Disp_mag,N_points); 

