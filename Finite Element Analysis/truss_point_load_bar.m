% ENME302
% Thomas Mitchell
% Truss point load frame

clear
clc
close all
format compact
format long

% ------------------------------ Parameters -------------------------------

Disp_mag = 100; % Magnification factor
N_points = 500; % Number of points on plot
Num_elements = 10; % Total number of elements in truss
E = 200*10^9; % Young's modulus in Pa
D1 = 0.150; % Bar outside diameter in m
D2 = 0.140; % Bar inside diameter in m
A = pi/4*(D1^2-D2^2); % Cross-sectional area of bar in m^2

% Global forcing terms in N
Q_nodal = [20000; 0; 0; 0; 40000; 0; 0; 0; 20000; 0]; 

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

AM = zeros(10, 4, Num_elements); % Assembly matrices
K = zeros(2,2,Num_elements); % Local stiffness terms
Khat = zeros(4,4,Num_elements); % Global stiffness terms
Lambda = zeros(2,4,Num_elements); % Transformation matrices
KG = zeros(10,10,Num_elements); % Structural contribution
D = zeros(4,1,Num_elements); % Global deflections
d = zeros(2,1,Num_elements); % Local deflections
f = zeros(2,1,Num_elements); % Local nodal forces
F = zeros(4,1,Num_elements); % Global nodal forces

% -------------------------- Assembly Matrices ----------------------------

AM(:,:,1) = [0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,2) = [1 0 0 0
             0 1 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,3) = [0 0 0 0
             0 0 0 0
             0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,4) = [1 0 0 0
             0 1 0 0
             0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,5) = [1 0 0 0
             0 1 0 0
             0 0 0 0
             0 0 0 0
             0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,6) = [0 0 0 0
             0 0 0 0
             0 0 1 0
             0 0 0 1
             1 0 0 0
             0 1 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0];

AM(:,:,7) = [0 0 0 0
             0 0 0 0
             1 0 0 0
             0 1 0 0
             0 0 0 0
             0 0 0 0
             0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0];

AM(:,:,8) = [0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             1 0 0 0
             0 1 0 0
             0 0 1 0
             0 0 0 1
             0 0 0 0
             0 0 0 0];

AM(:,:,9) = [0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             1 0 0 0
             0 1 0 0
             0 0 0 0
             0 0 0 0
             0 0 1 0
             0 0 0 1];

AM(:,:,10) = [0 0 0 0
              0 0 0 0
              0 0 0 0
              0 0 0 0
              0 0 0 0
              0 0 0 0
              0 0 1 0
              0 0 0 1
              1 0 0 0
              0 1 0 0];

% ---------------------------- Stiffness Terms ----------------------------

for i=1:Num_elements
    K(:,:,i) = local_bar(E,A,L(i)); % Element stiffness matrix in local coords
    [Khat(:,:,i), Lambda(:,:,i)] = global_bar(K(:,:,i),a(i)); % Global coords
    KG(:,:,i) = assembly_bar(Khat(:,:,i),AM(:,:,i)); % Structural contribution
end 

% --------------------------- Structural Terms ----------------------------

KG_total = KG(:,:,1);

for i=2:Num_elements
    KG_total = KG_total + KG(:,:,i); % Structural stiffness matrix
end 

q = KG_total\Q_nodal; % Structural displacements 

% ---------------------------- Displacements ------------------------------

for i=1:Num_elements
    D(:,i) = AM(:,:,i)'*q; % Global deflections
    d(:,i) = Lambda(:,:,i)*D(:,:,i); % Local deflections
    f(:,i) = K(:,:,i)*d(:,:,i); % Local nodal forces
    F(:,i) = Khat(:,:,i)*D(:,:,i); % Global nodal forces
    abs(F(1,i));
end 

% --------------------- Reactions and Deflections -------------------------

R1 = F(1:2,1,1); % Reaction forces at left support
R2 = F(3:4,1,2)+F(1:2,1,3); % Reaction forces at right support
D_tip = 1000*D(1:2,1,10); % Deflection components at top of truss

for i=1:Num_elements
    axial_stress = abs(f(1,1,i)/A); % Axial stress in element
    total_stress = axial_stress; % Absolute total stress
end

% --------------------------------- Plot ----------------------------------

hold on
for i=1:Num_elements
    Plot_deflected_bar(node1(i,1),node1(i,2),node2(i,1),node2(i,2),D(:,:,i),Disp_mag);
end

xlabel('Width (m)') 
ylabel('Height (m)')
title('Deflected vs Orignal shape for Bar Element Assumption')
legend({'Original Structure','Scaled Deflected Structure'},'Location','northeast')