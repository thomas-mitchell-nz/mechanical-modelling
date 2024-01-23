function [ Deflected_XG, Deflected_YG ] = Plot_deflected_shape ...
    (node1XG,node1YG,node2XG,node2YG,d_e,L,alpha,Disp_mag,N_points)
% PLOT_DEFLECTED_SHAPE Generate a plot for the deflected shape of a generic
% frame element

% Vector of points
x_e = linspace(0,L,N_points);

% Axial shape functions
psi1 = (1-x_e/L);
psi2 = x_e/L;

% Transverse shape functions
N1 = 1-3*x_e.^2/L^2+2*x_e.^3/L^3;
N2 = x_e.^3/L^2-2*x_e.^2/L+x_e;
N3 = 3*x_e.^2/L^2-2*x_e.^3/L^3;
N4 = x_e.^3/L^2-x_e.^2/L;

% Element axial displacement 
u = psi1*d_e(1) + psi2*d_e(4);

% Element transverse displacement
v = N1*d_e(2)+N2*d_e(3)+N3*d_e(5)+N4*d_e(6);

% Transfrom displacements from local to global coords
Deflections_XG = u*cosd(alpha)-v*sind(alpha);
Deflections_YG = u*sind(alpha)+v*cosd(alpha);

% Undeflected baseline of the element
Undeflected_baseline_XG = linspace(node1XG,node2XG,N_points);
Undeflected_baseline_YG = linspace(node1YG,node2YG,N_points);

% Deflected position of all points within the element
Deflected_XG = Undeflected_baseline_XG+Disp_mag*Deflections_XG;
Deflected_YG = Undeflected_baseline_YG+Disp_mag*Deflections_YG;

% Plot
plot(Undeflected_baseline_XG,Undeflected_baseline_YG,'b.-')
hold on
plot(Deflected_XG,Deflected_YG,'r.-')

end