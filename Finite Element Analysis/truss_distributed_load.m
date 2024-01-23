% ENME302
% Thomas Mitchell
% Truss distributed load

clear
clc
close all
format compact

% ------------------------------ Parameters -------------------------------

Disp_mag = 100; % Magnification factor
N_points = 10; % Number of points on plot
Num_elements = 10; % Total number of elements in truss
E = 200*10^9; % Young's modulus in Pa
D1 = 0.150; % Frame outside diameter in m
D2 = 0.140; % Frame inside diameter in m
A = pi/4*(D1^2-D2^2); % Cross-sectional area of frame in m^2
I = 1/4*pi*((D1/2)^4-(D2/2)^4);  % 2nd moment of area of frame
c = 75*10^-3; % Distance to the outermost fibre for bending calculations 

% Wind loads in Pa:
low_pressure = 0.61 * 10^3;
medium_pressure = 0.82*10^3;
high_pressure = 1.16*10^3;
very_high_pressure = 1.50*10^3;

% Wind loads in Pa/m
frame_depth = 10; % m
wbar_low = low_pressure * frame_depth;
wbar_medium = medium_pressure * frame_depth;
wbar_high = high_pressure * frame_depth;
wbar_very_high = very_high_pressure * frame_depth;

% Yield properties:
yield_stress = 280*10^6; % Pa
FoS = 2; % factor of sagety
max_allowable_stress = yield_stress/FoS;
max_wbar = wbar_low*3.53; % Guess and check
max_pressure = max_wbar / frame_depth; % Pa
max_wind_speed_ms = sqrt(max_pressure/0.6); % m/s
max_wind_speed_kmh = max_wind_speed_ms * 3.6;

active_wbar = wbar_low; % current wbar in use by program

% Element lengths in m:
L_h = 2;
L_v = 2.5;
L_d = sqrt(L_h^2+L_v^2);
L = [L_v L_d L_v L_h L_v L_d L_v L_h L_v L_d];

% Transformations angles in degrees:
a_h = 0;
a_v = 90;
a_d = -atand(L_v/L_h);
a = [a_v a_d a_v a_h a_v a_d a_v a_h a_v a_d];

% Global position of node 1 (x,y) for each element:
node1 = [0 0;
         0 L_v;
         L_h 0;
         0 L_v;
         0 L_v;
         0 2*L_v;
         L_h L_v;
         0 2*L_v;
         0 2*L_v;
         0 3*L_v];

% Global position of node 2 (x,y) for each element:
node2 = [0 L_v;
         L_h 0;
         L_h L_v;
         L_h L_v;
         0 2*L_v;
         L_h L_v;
         L_h 2*L_v;
         L_h 2*L_v;
         0 3*L_v;
         L_h 2*L_v];

% ---------------------------- Preallocations -----------------------------

AM = zeros(15, 6, Num_elements); % Assembly matrices
K = zeros(6,6,Num_elements); % Local stiffness terms
Khat = zeros(6,6,Num_elements); % Global stiffness terms
Lambda = zeros(6,6,Num_elements); % Transformation matrices
KG = zeros(15,15,Num_elements); % Structural contribution
D = zeros(6,1,Num_elements); % Global deflections
d = zeros(6,1,Num_elements); % Local deflections
f = zeros(6,1,Num_elements); % Local nodal forces
F = zeros(6,1,Num_elements); % Global nodal forces

% -------------------------- Assembly Matrices ----------------------------

AM(:,:,1) = [0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,2) = [1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,3) = [0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,4) = [1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,5) = [1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,6) = [0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1;
             1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,7) = [0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 1 0 0;
             0 0 0 0 1 0; 
             0 0 0 0 0 1;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,8) = [0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 1 0 0;
             0 0 0 0 1 0; 
             0 0 0 0 0 1;
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0];

AM(:,:,9) = [0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             1 0 0 0 0 0; 
             0 1 0 0 0 0; 
             0 0 1 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 1 0 0; 
             0 0 0 0 1 0; 
             0 0 0 0 0 1];

AM(:,:,10) = [0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 0 0 0;
              0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 0 0 0; 
              0 0 0 1 0 0;
              0 0 0 0 1 0; 
              0 0 0 0 0 1;
              1 0 0 0 0 0; 
              0 1 0 0 0 0; 
              0 0 1 0 0 0];

% ---------------------------- Stiffness Terms ----------------------------

for i=1:Num_elements
    K(:,:,i) = local_frame(E,I,A,L(i)); % Element stiffness matrix in local coords
    [Khat(:,:,i), Lambda(:,:,i)] = global_frame(K(:,:,i),a(i)); % Global coords
    KG(:,:,i) = assembly_frame(Khat(:,:,i),AM(:,:,i)); % Structural contribution
end 

% --------------------------- Equivalent Loads ----------------------------

% Uniformly distributed load:
[f5_UDL,F5_UDL,Q5_UDL] = equiv_uniform_load(-active_wbar,L(5),Lambda(:,:,5),AM(:,:,5));
[f9_UDL,F9_UDL,Q9_UDL] = equiv_uniform_load(-active_wbar,L(9),Lambda(:,:,9),AM(:,:,9));

% --------------------------- Structural Terms ----------------------------

KG_total = KG(:,:,1);
Q_total = Q5_UDL + Q9_UDL;

for i=2:Num_elements
    KG_total = KG_total + KG(:,:,i); % Structural stiffness matrix
end 
q = KG_total\Q_total; % Structural displacements 

% ---------------------------- Displacements ------------------------------

for i=1:Num_elements
    D(:,i) = AM(:,:,i)'*q; % Global deflections
    d(:,i) = Lambda(:,:,i)*D(:,:,i); % Local deflections
    f(:,i) = K(:,:,i)*d(:,:,i); % Local nodal forces
    F(:,i) = Khat(:,:,i)*D(:,:,i); % Global nodal forces
end 

% --------------------- Reactions and Deflections -------------------------

R1=F(1:3,1,1); % Reaction forces at left support
R2=F(4:6,1,2)+F(1:3,1,3); % Reaction forces at right support

M1=-20000*2.5-40000*2*2.5-20000*3*2.5+R2(2)*2+R2(3); % Moment at left support
M2=-20000*2.5-40000*2*2.5-20000*3*2.5-R1(2)*2+R1(3); % Moment at right support

D_tip = D(1:3,1,10); % Deflection components at top of truss

for i=1:Num_elements
    axial_stress = abs(f(1,1,i)/A); % Axial stress in element
    bending_stress = (max(abs(f(3,1,i)),abs(f(6,1,i)))*c/I); % Max bending stress in element
    total_stress = axial_stress + bending_stress; % Absolute total stress
end

% --------------------------------- Plot ----------------------------------

hold on
for i=1:Num_elements
    Plot_deflected_shape(node1(i,1),node1(i,2),node2(i,1),node2(i,2),d(:,:,i),L(i),a(i),Disp_mag,N_points);
end

xlabel('Width (m)') 
ylabel('Height (m)')
title('Deflected vs Orignal shape for Max Distributed Wind Loading')
legend({'Original Structure','Scaled Deflected Structure'},'Location','northeast')